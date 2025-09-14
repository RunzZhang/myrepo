import uproot
import pandas as pd
import os, time


class ReadRoot:
    def __init__(self):
        self.base_path = "/data/runzezhang/result/TN_sims_D/"
        self.plot_path = '/data/runzezhang/result/TN_sims_D/plot/'

        # self.base_path = "/data/runzezhang/result/TN_box/"
        # self.plot_path = '/data/runzezhang/result/TN_box/plot/'

        # self.false_1 = "AmLi_1E7_false1.csv"
        # self.false_2 = "AmLi_1E7_false2.csv"
        # self.signal = "AmLi_1E7_sig.csv"
        # self.false_1_mid = "AmLi_1E7_false1_mid.csv"
        # self.false_2_mid = "AmLi_1E7_false2_mid.csv"
        # self.signal_mid = "AmLi_1E7_sig_mid.csv"
        #
        # self.false_1_path = self.base_path + self.false_1
        # self.false_2_path = self.base_path + self.false_2
        # self.false_1_path_mid = self.base_path + self.false_1_mid
        # self.false_2_path_mid = self.base_path + self.false_2_mid
        # self.signal_path_mid = self.base_path + self.signal_mid
        # self.signal_path = self.base_path + self.signal

        self.filepath = self.base_path + "dmx_AmLi_1E7.root"
        # self.filepath = self.base_path + "dmx_AmLi.root" # test
        self.tree_name = "tree"  # Assuming your TTree is named "tree"

        # Define the columns you want to read and write
        # self.selected_columns = ["Event", "name", "Parent ID", "Track ID", "Step ID", "X/mm", "PreKinetic/MeV","PostKinetic/MeV"
        #                          "Recoiled/MeV", "Volume", "Process"]
        self.selected_columns = ["Event", "name", "Parent ID", "Track ID", "Step ID", "X/mm", "PreKinetic/MeV",
                                 "PostKinetic/MeV",
                                 "Recoiled/MeV", "Volume", "Process"]

    def chunk_and_write_root(self, num_chunks=20, output_dir=None):
        if output_dir is None:
            output_dir = os.path.join(self.base_path, "chunked_root_files_AmLi_SN")
        os.makedirs(output_dir, exist_ok=True)

        with uproot.open(self.filepath) as file:
            tree = file[self.tree_name]
            total_entries = int(tree.num_entries)
            print(f"Total entries in original file: {total_entries}")

            # Calculate entries per chunk
            entries_per_chunk = total_entries // num_chunks
            # Ensure all entries are covered, handle remainder
            entry_steps = [entries_per_chunk] * (num_chunks - 1)
            entry_steps.append(total_entries - sum(entry_steps))

            print(f"Splitting into {num_chunks} chunks.")
            print(f"Entries per chunk (approx): {entries_per_chunk}")

            # Get the schema (data types) of the selected columns for writing
            # You need to open the tree to get the types. We can get it from a small head.
            # Convert to dictionary format as required by uproot.recreate
            branch_types = {col: tree[col].interpretation.numpy_dtype for col in self.selected_columns}

            chunk_num = 0
            start_entry = 0
            # Iterate through the ROOT file in chunks
            for arrays in tree.iterate(expressions=self.selected_columns, library="pd", entry_start=0,
                                       entry_stop=total_entries, step_size=entries_per_chunk):
                arrays['name'] = arrays['name'].astype(str)
                arrays['Volume'] = arrays['Volume'].astype(str)
                arrays['Process'] = arrays['Process'].astype(str)

                chunk_num += 1
                output_filename = os.path.join(output_dir, f"dmx_AmLi_1E7_2_{chunk_num}.root")

                print(f"Processing chunk {chunk_num} (entries {start_entry} to {start_entry + len(arrays) - 1})")

                # Convert the pandas DataFrame to a dictionary of NumPy arrays
                # This is the format uproot.recreate expects for TTree branches
                data_to_write = {col: arrays[col].values for col in self.selected_columns}

                try:
                    with uproot.recreate(output_filename, compression=uproot.LZ4(level=4)) as output_file:
                        output_file[self.tree_name] = data_to_write
                    print(f"Successfully wrote {len(arrays)} entries to {output_filename}")
                except Exception as e:
                    print(f"Error writing chunk {chunk_num} to {output_filename}: {e}")

                start_entry += len(arrays)
                time.sleep(5)

        print(f"\nFinished chunking the ROOT file into {chunk_num} files in {output_dir}.")


# Example usage:
if __name__ == "__main__":
    reader = ReadRoot()
    reader.chunk_and_write_root(num_chunks=100)