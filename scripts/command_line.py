import sys
import os
from os.path import basename, isfile
from .seeker import SeekerFasta

# Constants ------------------------------------------------------------------------------------------------------------
MGM_USAGE = "Usage:\npredict-metagenome <input fasta file>"


# Terminal functions ---------------------------------------------------------------------------------------------------
def predict_metagenome():
    if len(sys.argv) != 2:
        print(MGM_USAGE)
        sys.exit(1)

    path = sys.argv[1]
    if not isfile(path):
        print("{} is not a valid file path.\n{}".format(path, MGM_USAGE))
        sys.exit(1)

    print("Reading Fasta at {}...".format(path), file=sys.stderr)
    seeker_fasta = SeekerFasta(path, load_seqs=False)
    print("name\tprediction\tscore")
    
    predictions = seeker_fasta.phage_or_bacteria()
    abs_path= os.path.dirname(os.path.abspath(sys.argv[1]))
    file_base=(os.path.split(sys.argv[1])[1]).strip(".fa")
    output=abs_path + "/" + "seeker_" + file_base + "_output.txt"
    with open(output, "a+") as f:
        for i in predictions:
            f.write(i + "\n")
    print("\n".join(seeker_fasta.phage_or_bacteria()))
