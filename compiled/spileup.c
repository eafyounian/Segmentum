
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>

#undef assert
#define assert(cond) if (!(cond)) { fprintf(stderr, "Assert '%s' failed in %s, line %d.\nError encountered while parsing line %d of pileup.\n", #cond, __FILE__, __LINE__, line_num); exit(-1); } 

#define MAX_ALLELE_LEN 32
#define MAX_ALLELES 100
#define MAX_SAMPLES 5000
#define MAX_TABS_ON_LINE (5 + 4*MAX_SAMPLES)

typedef struct {
	char seq[MAX_ALLELE_LEN];
	int high_count;			// High quality reads
	int low_count;			// Low quality reads
} Allele;

typedef struct {
	Allele alleles[MAX_ALLELES];
	int total;
} AlleleList;

int min_mapq = 0;
uint64_t line_num = 0;       // Line number

void count_allele(char* allele, int len, char quality, AlleleList* list) {
	quality = quality - 33;

	// Zero terminate temporarily.
	char ochar = allele[len]; allele[len] = 0;
	
	for (int i = 0; i < list->total; i++) {
		if (strcmp(allele, list->alleles[i].seq) == 0) {
			if (quality >= min_mapq)
				list->alleles[i].high_count += 1;
			else
				list->alleles[i].low_count += 1;
			goto cleanup;
		}
	}
	
	// Haven't seen this allele yet. Add it to the list unless there's
	// already too many alleles at this position.
	if (list->total >= MAX_ALLELES) goto cleanup;
	list->alleles[list->total].high_count = (quality >= min_mapq);
	list->alleles[list->total].low_count = (quality < min_mapq);
	strncpy(list->alleles[list->total].seq, allele, MAX_ALLELE_LEN); 
	list->total += 1;
	
cleanup:
	allele[len] = ochar;   // Undo zero termination
}

void print_alleles(AlleleList* list) {
	if (list->total == 0) return;
	int i;
	for (i = 0; i < list->total - 1; i++) {
		printf("%s %d %d ", list->alleles[i].seq, list->alleles[i].high_count,
			list->alleles[i].low_count);
	}
	assert(i < MAX_ALLELES);
	printf("%s %d %d", list->alleles[i].seq, list->alleles[i].high_count,
		list->alleles[i].low_count);
}

int parse_pileup(char* bases, char* quality, AlleleList* alleles) {
	// First we need to convert everything into lowercase.
	for (int i = 0; bases[i]; i++) {
		if (bases[i] == ',') {
			bases[i] = '.';
		} else {
			bases[i] = toupper(bases[i]);
		}
	}
	
	int j = 0;
	for (int i = 0; bases[i]; i++) {
		if (bases[i] == '^') {
			i++;
		} else if (bases[i] == '*' || bases[i] == '$' || bases[i] == '<' ||
			bases[i] == '>') {
			// Ignore.
		} else if (bases[i+1] == '+' || bases[i+1] == '-') {
			char* allele = 0;  // First char after number stored here.
			int len = strtol(bases + i + 2, &allele, 10);
			if (len < MAX_ALLELE_LEN) {
				allele[-1] = bases[i+1]; allele[-2] = bases[i];
				allele -= 2; len += 2;  // Include pre-base and sign
				count_allele(allele, len, quality[j++], alleles);
			}
			i += len;
		} else if (bases[i] == 'N') {
			// Unknown bases should not count towards the total.
		} else {
			//assert(bases[i] == 'A' || bases[i] == 'C' || bases[i] == 'G' || bases[i] == 'T' || bases[i] == '.');
			count_allele(bases + i, 1, quality[j++], alleles);
		}
	}
}

int main(int argc, char** argv) {
	if (argc != 3) {
		printf("Usage: spileup <min_alt_reads> <min_mapq>\n");
		return -1;
	}

	// Parse command line arguments.
	int min_alt_reads = strtol(argv[1], NULL, 10);
	if (errno == EINVAL || min_alt_reads < 0) {
		fprintf(stderr, "Invalid min_alt_reads argument.\n");
		return -1;
	}
	fprintf(stderr, "Pre-filtering variants with less than %d alt reads...\n",
		min_alt_reads);

	min_mapq = strtol(argv[2], NULL, 10);
	if (errno == EINVAL || min_mapq < 0) {
		fprintf(stderr, "Invalid min_mapq argument.\n");
		return -1;
	}
	fprintf(stderr, "Ignoring reads with MAPQ below %d...\n", min_mapq);
	
	size_t line_buf_size = 1024 * 1024;
	char* line = malloc(line_buf_size);
	
	int* tabs = malloc(sizeof(int) * MAX_TABS_ON_LINE);  // Tab locations
	AlleleList* sample_alleles = malloc(sizeof(AlleleList) * MAX_SAMPLES);
	
	// Expected format is:
	// CHROMOSOME, POSITION, REF, [READ_COUNT, READ_BASES, BASEQ, MAPQ]+
	while (getline(&line, &line_buf_size, stdin) != -1) {
		line_num += 1;

		// Find all tab characters.
		int ntabs = 0;
		for (int i = 0; line[i]; i++) {
			if (line[i] == '\t') {
				assert(ntabs < MAX_TABS_ON_LINE);
				tabs[ntabs++] = i; line[i] = '\0';
			}
		}
		if (ntabs < 3) continue;

		int S = 0, t = 2;
		while (t < ntabs) {
			assert(S < MAX_SAMPLES);
			sample_alleles[S].total = 0;
			assert(isdigit(line[tabs[t]+1]));

			// Note: Before samtools-1.1, if there were no reads overlapping
			// a base, samtools mpileup would output 3 columns "0\t*\t*"
			// instead of four. In samtools-1.1, output is always 4 columns.
			parse_pileup(line + tabs[t+1] + 1, line + tabs[t+3] + 1,
				&sample_alleles[S]);
			t += 4; S += 1;
		}
		
		// If all samples only show reference alleles, don't print anything.
		int best_alt_reads = 0;
		for (int s = 0; s < S; s++) {
			AlleleList* alleles = &sample_alleles[s];
			int alt_alleles = 0;
			for (int a = 0; a < alleles->total; a++) {
				assert(a < MAX_ALLELES);
				if (alleles->alleles[a].seq[0] == '.' &&
					alleles->alleles[a].seq[1] == '\0') continue;
				if (alleles->alleles[a].high_count > best_alt_reads)
					best_alt_reads = alleles->alleles[a].high_count;
			}
		}
		
		if (best_alt_reads >= min_alt_reads) {
			printf("%s\t%s\t%s", line, line + tabs[0]+1, line + tabs[1]+1); 
			for (int s = 0; s < S; s++) {
				printf("\t");
				print_alleles(&sample_alleles[s]);
			}
			printf("\n");
			fflush(stdout);
		}
	}
	return 0;
}


