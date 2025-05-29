#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_LEN 2048
#define SUBSTR_LEN 1024

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s [MEM_FILE]\n", argv[0]);
        fprintf(stderr, "[MEM_FILE] is the output of running ropebwt3 with KeBaB fragments\n");
        return 1;
    }

    FILE *fp = fopen(argv[1], "r");
    if (!fp) {
        perror("Error opening file");
        return 1;
    }

    char line[LINE_LEN];
    while (fgets(line, sizeof(line), fp)) {
        char seq[SUBSTR_LEN];
        int start, end, mem_start, mem_end, occ;
        if (sscanf(line, "%[^:]:%d-%d\t%d\t%d\t%d\n", seq, &start, &end, &mem_start, &mem_end, &occ) == 6) {
            int offset = start - 1;
            printf("%s\t%d\t%d\t%d\n", seq, mem_start + offset, mem_end + offset, occ);
        }
    }

    fclose(fp);
    return 0;
}