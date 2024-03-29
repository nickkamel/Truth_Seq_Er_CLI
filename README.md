# Truth-Seq-Er

Truth-Seq-Er is a multi-objective evolutionary algorithm that designs ribozyme-based logic gates (ribogates) implementing target 1, 2, and 3-input functions.

## Requirements
- Python 3.x
- Numpy

## Setup

- Clone this repository.
- Place the four folding binaries (**RNAfold.exe**, **RNAcofold.exe**, **RNAfold300.exe**, **RNAcofold300.exe**) in the same directory as the main script (truth_seq_er.py). These binaries were generated by compiling the ViennaRNA package with the **MAXLOOP** parameter changed from 30 to 300, which is necessary for correct structure prediction of larger ribogates. Note that slight changes have been made to these binaries to print out the base-pairing probability matrix to the console instead of generating .ps files.
- (Optional) To use the ViennaRNA Python bindings instead of the provided folding binaries, modify the vienna.py script to interface with the bindings.

## Usage

Run the script with the following command:

`python truth_seq_er.py --num_inputs=[NUM INPUTS] --truth_vector=[TRUTH VECTOR]`

For example:

`python truth_seq_er.py --num_inputs=3 --truth_vector=10000001`

Here, the truth vector corresponds to the last column of the truth table for the specified function.

Truth-Seq-Er will generate a diverse population of ribogates implementing the target function. The designs are stored in a CSV file named **ribogate_designs_[TARGET FUNCTION NAME].csv**. Each row in the file corresponds to a designed ribogate and contains the sequence of the ribogate and input strands. A second CSV file named **ribogate_structures_[TARGET FUNCTION NAME].csv** provides the predicted secondary structure in dot-bracket format for each state of the ribogates.

## Note
Ribogate design is computationally intensive; designing a 3-input gate may take several hours on a multi-core machine.
