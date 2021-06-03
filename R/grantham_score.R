# Return Grantham score for protein changes (extracted from annotations)

grantham_score <- function(annotation)
{
    # Function based on: https://github.com/ashutoshkpandey/Annotation/blob/master/Grantham_score_calculator.py

    # Granthan score matrix (use upper diagonal only and transpose)
    #         Arg Leu Pro Thr Ala Val Gly Ile Phe Tyr Cys His Gln Asn Lys Asp Glu Met Trp     
    #         R   L   P   T   A   V   G   I   F   Y   C   H   Q   N   K   D   E   M   W       
    # Ser S   110 145 74  58  99  124 56  142 155 144 112 89  68  46  121 65  80  135 177 S   Ser
    # Arg R   0   102 103 71  112 96  125 97  97  77  180 29  43  86  26  96  54  91  101 R   Arg
    # Leu L   0   0   98  92  96  32  138 5   22  36  198 99  113 153 107 172 138 15  61  L   Leu
    # Pro P   0   0   0   38  27  68  42  95  114 110 169 77  76  91  103 108 93  87  147 P   Pro
    # Thr T   0   0   0   0   58  69  59  89  103 92  149 47  42  65  78  85  65  81  128 T   Thr
    # Ala A   0   0   0   0   0   64  60  94  113 112 195 86  91  111 106 126 107 84  148 A   Ala
    # Val V   0   0   0   0   0   0   109 29  50  55  192 84  96  133 97  152 121 21  88  V   Val
    # Gly G   0   0   0   0   0   0   0   135 153 147 159 98  87  80  127 94  98  127 184 G   Gly
    # Ile I   0   0   0   0   0   0   0   0   21  33  198 94  109 149 102 168 134 10  61  I   Ile
    # Phe F   0   0   0   0   0   0   0   0   0   22  205 100 116 158 102 177 140 28  40  F   Phe
    # Tyr Y   0   0   0   0   0   0   0   0   0   0   194 83  99  143 85  160 122 36  37  Y   Tyr
    # Cys C   0   0   0   0   0   0   0   0   0   0   0   174 154 139 202 154 170 196 215 C   Cys
    # His H   0   0   0   0   0   0   0   0   0   0   0   0   24  68  32  81  40  87  115 H   His
    # Gln Q   0   0   0   0   0   0   0   0   0   0   0   0   0   46  53  61  29  101 130 Q   Gln
    # Asn N   0   0   0   0   0   0   0   0   0   0   0   0   0   0   94  23  42  142 174 N   Asn
    # Lys K   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   101 56  95  110 K   Lys
    # Asp D   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   45  160 181 D   Asp
    # Glu E   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   126 152 E   Glu
    # Met M   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   67  M   Met

    Grantham_matrix <- matrix(c(  0, 110, 145,  74,  58,  99, 124,  56, 142, 155, 144, 112,  89,  68,  46, 121,  65,  80, 135, 177,
                                110,   0, 102, 103,  71, 112,  96, 125,  97,  97,  77, 180,  29,  43,  86,  26,  96,  54,  91, 101,
                                145, 102,   0,  98,  92,  96,  32, 138,   5,  22,  36, 198,  99, 113, 153, 107, 172, 138,  15,  61,
                                 74, 103,  98,   0,  38,  27,  68,  42,  95, 114, 110, 169,  77,  76,  91, 103, 108,  93,  87, 147,
                                 58,  71,  92,  38,   0,  58,  69,  59,  89, 103,  92, 149,  47,  42,  65,  78,  85,  65,  81, 128,
                                 99, 112,  96,  27,  58,   0,  64,  60,  94, 113, 112, 195,  86,  91, 111, 106, 126, 107,  84, 148,
                                124,  96,  32,  68,  69,  64,   0, 109,  29,  50,  55, 192,  84,  96, 133,  97, 152, 121,  21,  88,
                                 56, 125, 138,  42,  59,  60, 109,   0, 135, 153, 147, 159,  98,  87,  80, 127,  94,  98, 127, 184,
                                142,  97,   5,  95,  89,  94,  29, 135,   0,  21,  33, 198,  94, 109, 149, 102, 168, 134,  10,  61,
                                155,  97,  22, 114, 103, 113,  50, 153,  21,   0,  22, 205, 100, 116, 158, 102, 177, 140,  28,  40,
                                144,  77,  36, 110,  92, 112,  55, 147,  33,  22,   0, 194,  83,  99, 143,  85, 160, 122,  36,  37,
                                112, 180, 198, 169, 149, 195, 192, 159, 198, 205, 194,   0, 174, 154, 139, 202, 154, 170, 196, 215,
                                 89,  29,  99,  77,  47,  86,  84,  98,  94, 100,  83, 174,   0,  24,  68,  32,  81,  40,  87, 115,
                                 68,  43, 113,  76,  42,  91,  96,  87, 109, 116,  99, 154,  24,   0,  46,  53,  61,  29, 101, 130,
                                 46,  86, 153,  91,  65, 111, 133,  80, 149, 158, 143, 139,  68,  46,   0,  94,  23,  42, 142, 174,
                                121,  26, 107, 103,  78, 106,  97, 127, 102, 102,  85, 202,  32,  53,  94,   0, 101,  56,  95, 110,
                                 65,  96, 172, 108,  85, 126, 152,  94, 168, 177, 160, 154,  81,  61,  23, 101,   0,  45, 160, 181,
                                 80,  54, 138,  93,  65, 107, 121,  98, 134, 140, 122, 170,  40,  29,  42,  56,  45,   0, 126, 152,
                                135,  91,  15,  87,  81,  84,  21, 127,  10,  28,  36, 196,  87, 101, 142,  95, 160, 126,   0,  67,
                                177, 101,  61, 147, 128, 148,  88, 184,  61,  40,  37, 215, 115, 130, 174, 110, 181, 152,  67,   0), 
                                nrow=20, byrow=TRUE)
        rownames(Grantham_matrix) <- c("S", "R", "L", "P", "T", "A", "V", "G", "I", "F", "Y", "C", "H", "Q", "N", "K", "D", "E", "M", "W")
        colnames(Grantham_matrix) <- rownames(Grantham_matrix)

    # Get amino acid change part of annotation
    last_strsplit <- function(x, split) {
        x_split <- strsplit(x, split)[[1]]
        if(length(x_split) == 0) {
            return(".")
        }
        return(as.vector(unlist(x_split[length(x_split)])))
    }
    annotation <- sapply(annotation, function(x){last_strsplit(x, ":")})
    
    # Remove annotations with no protein changes
    annotation[!startsWith(annotation, "p.")] <- "."
    
    # Remove "p." and any commas
    annotation <- gsub("p.", "", annotation)
    annotation <- gsub(",", "", annotation)

    # Replace entries with deletions, frameshipts, unknown or wholegene by "." as will not return Grantham score
    annotation[endsWith(annotation, "del")] <- "."
    annotation[endsWith(annotation, "fs")] <- "."

    # Translate from three letter amino acid codes to one letter
    annotation <- gsub("Ala", "A", annotation)
    annotation <- gsub("Arg", "R", annotation)
    annotation <- gsub("Asn", "N", annotation)
    annotation <- gsub("Asp", "D", annotation)
    annotation <- gsub("Cys", "C", annotation)
    annotation <- gsub("Gln", "Q", annotation)
    annotation <- gsub("Glu", "E", annotation)
    annotation <- gsub("Gly", "G", annotation)
    annotation <- gsub("His", "H", annotation)
    annotation <- gsub("Ile", "I", annotation)
    annotation <- gsub("Leu", "L", annotation)
    annotation <- gsub("Lys", "K", annotation)
    annotation <- gsub("Met", "M", annotation)
    annotation <- gsub("Phe", "F", annotation)
    annotation <- gsub("Pro", "P", annotation)
    annotation <- gsub("Pyl", "O", annotation)
    annotation <- gsub("Ser", "S", annotation)
    annotation <- gsub("Sec", "U", annotation)
    annotation <- gsub("Thr", "T", annotation)
    annotation <- gsub("Trp", "W", annotation)
    annotation <- gsub("Tyr", "Y", annotation)
    annotation <- gsub("Val", "V", annotation)
    annotation <- gsub("Asx", "B", annotation)
    annotation <- gsub("Glx", "Z", annotation)
    annotation <- gsub("Xaa", "X", annotation)
    annotation <- gsub("Xle", "J", annotation)

    # "X", "B", "Z", "J" represent ambiguous amino acids, so remove these
    annotation[startsWith(annotation, "X") | endsWith(annotation, "X")] <- "."
    annotation[startsWith(annotation, "B") | endsWith(annotation, "B")] <- "."
    annotation[startsWith(annotation, "Z") | endsWith(annotation, "Z")] <- "."
    annotation[startsWith(annotation, "J") | endsWith(annotation, "J")] <- "."

    # Get before and after proteins and calculate Grantham score
    protein_before <- substr(annotation, 1, 1)
    protein_after <- sapply(annotation, function(x){substr(x, nchar(x), nchar(x))})
    Grantham_scores <- rep(NA, length(annotation))
    for (ii in 1:length(annotation)) {
        if (protein_before[ii] %in% rownames(Grantham_matrix) & protein_after[ii] %in% colnames(Grantham_matrix)) {
            Grantham_scores[ii] <- Grantham_matrix[protein_before[ii], protein_after[ii]]
        }
    }
    return(Grantham_scores)
}


