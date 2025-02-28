#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#define MAX_CPG 1000
#define MAX_READS 1000

typedef struct {
    char chrom[100];
    int pos;
} CpGSite;

CpGSite cpg_sites[MAX_CPG];
int cpg_count = 0;

// Matrix to store likelihoods
// Rows = reads, columns = CpG sites
int likelihood_matrix[MAX_READS][MAX_CPG];
char read_names[MAX_READS][100];
int read_count = 0;

// Function to load CpG sites from file
void load_cpg_sites(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening CpG sites file\n");
        exit(1);
    }
    while (fscanf(file, "%99[^:]:%d\n", cpg_sites[cpg_count].chrom, &cpg_sites[cpg_count].pos) == 2) {
        cpg_count++;
    }
    fclose(file);
    printf("Loaded %d CpG sites\n", cpg_count);
}

// Function to find index of CpG site in the list
int find_cpg_index(const char *chrom, int pos) {
    for (int i = 0; i < cpg_count; i++) {
        if (strcmp(cpg_sites[i].chrom, chrom) == 0 && cpg_sites[i].pos == pos) {
            return i;
        }
    }
    return -1;
}

// Function to extract and store methylation likelihoods, ensuring indels are correctly handled
void extract_methylation(bam1_t *aln, bam_hdr_t *header) {
    const char *read_name = bam_get_qname(aln);
    strcpy(read_names[read_count], read_name);

    // Get MM and ML tags
    uint8_t *mm = bam_aux_get(aln, "MM");
    uint8_t *ml = bam_aux_get(aln, "ML");
    if (!mm || !ml) return;
    mm++; // Skip type character
    ml++; // Skip type character

    // CIGAR string
    uint32_t *cigar = bam_get_cigar(aln);
    int ref_pos = aln->core.pos;
    int read_pos = 0;
    int cpg_index;

    // Parse CIGAR to map read to reference positions
    for (int k = 0; k < aln->core.n_cigar; k++) {
        int op_len = bam_cigar_oplen(cigar[k]);
        int op = bam_cigar_op(cigar[k]);

        switch (op) {
            case BAM_CMATCH: // Alignment match (can be mismatch)
            case BAM_CEQUAL: // Sequence match
            case BAM_CDIFF:  // Sequence mismatch
                for (int i = 0; i < op_len; i++) {
                    // Check if reference position is a target CpG
                    cpg_index = find_cpg_index(header->target_name[aln->core.tid], ref_pos);
                    if (cpg_index >= 0) {
                        likelihood_matrix[read_count][cpg_index] = ml[read_pos];
                    }
                    ref_pos++;
                    read_pos++;
                }
                break;
            
            case BAM_CINS:  // Insertion to the reference
                read_pos += op_len;
                break;
            
            case BAM_CDEL:  // Deletion from the reference
                ref_pos += op_len;
                break;
            
            case BAM_CSOFT_CLIP:  // Soft clipping (read bases present)
                read_pos += op_len;
                break;
            
            case BAM_CHARD_CLIP:  // Hard clipping (read bases not present)
                break;
            
            case BAM_CREF_SKIP:  // Skipped region from the reference
                ref_pos += op_len;
                break;
            
            default:
                fprintf(stderr, "Unexpected CIGAR operation: %d\n", op);
                break;
        }
    }

    read_count++;
}

int main() {
    load_cpg_sites("/ifs/data/research/projects/juliet/tools/create_cpg_matrix/test/cpg_sites.txt");

    const char *filename = "/ifs/data/research/projects/juliet/tools/create_cpg_matrix/test/P50-A5.haplotagged.bam";
    samFile *infile = sam_open(filename, "r");
    if (infile == NULL) {
        fprintf(stderr, "Error opening BAM file\n");
        return 1;
    }

    bam_hdr_t *header = sam_hdr_read(infile);
    if (header == NULL) {
        fprintf(stderr, "Error reading header\n");
        sam_close(infile);
        return 1;
    }

    hts_idx_t *idx = sam_index_load(infile, filename);
    if (idx == NULL) {
        fprintf(stderr, "Error loading index\n");
        bam_hdr_destroy(header);
        sam_close(infile);
        return 1;
    }

    bam1_t *aln = bam_init1();

    // Query each CpG site
    for (int i = 0; i < cpg_count; i++) {
        char region[200];
        sprintf(region, "%s:%d-%d", cpg_sites[i].chrom, cpg_sites[i].pos, cpg_sites[i].pos);
        hts_itr_t *iter = sam_itr_querys(idx, header, region);

        if (iter == NULL) {
            fprintf(stderr, "No reads found in region: %s\n", region);
            continue;
        }

        while (sam_itr_next(infile, iter, aln) >= 0) {
            extract_methylation(aln, header);
        }

        hts_itr_destroy(iter);
    }

    // Output matrix
    printf("\nMethylation likelihood matrix:\n");
    printf("read_name\t");
    for (int j = 0; j < cpg_count; j++) {
        printf("%s_%d\t", cpg_sites[j].chrom, cpg_sites[j].pos);
    }
    printf("\n");
    for (int i = 0; i < read_count; i++) {
        printf("%s\t", read_names[i]);
        for (int j = 0; j < cpg_count; j++) {
            printf("%d\t", likelihood_matrix[i][j]);
        }
        printf("\n");
    }

    bam_destroy1(aln);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(infile);

    return 0;
}
