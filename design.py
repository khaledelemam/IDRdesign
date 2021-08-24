from designAlgo import design_sequences
import sys


def main():
    arguments = len(sys.argv) -1 
    if arguments == 2:  
        target_seq = sys.argv[1]
        seqs_number = sys.argv[2]
        query_seq = ''

    elif arguments == 3:
        target_seq =sys.argv[1]
        seqs_number = sys.argv[2]
        query_seq = sys.argv[3]

    
    else:
        print("Incorrect number of arguments")
        print("Please provide 1) target sequence 2) number of desired designed sequences 3) [Optional] query sequence")
        sys.exit(1)


    designed = design_sequences(target_seq, int(seqs_number), query_seq)
    # print(designed)


if __name__ == '__main__':
    main()