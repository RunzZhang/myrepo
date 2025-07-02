import uproot
import awkward as ak
base_address = "/data/runzezhang/result/TN_sims_D/"
root_name = "dmx_Cf_1E7.root"
input_file = base_address+root_name

tree_name = "tree"  # Replace with your tree name
columns = None  # or list of specific branches if you want to filter
num_parts = 10

# Open the full tree
with uproot.open(f"{input_file}:{tree_name}") as tree:
    total_entries = tree.num_entries
    chunk_size = total_entries // num_parts

    print(f"Total entries: {total_entries}")
    print(f"Each chunk: {chunk_size} entries")

    for i in range(num_parts):
        start = i * chunk_size
        end = (i + 1) * chunk_size if i < num_parts - 1 else total_entries

        print(f"Processing part {i}: entries {start} to {end}")

        # Read chunk
        array_chunk = tree.arrays(
            columns,
            entry_start=start,
            entry_stop=end,
            library="ak"
        )

        # Save to new root file
        output_file = f"dmx_Cf_1e7_part{i}.root"
        with uproot.recreate(base_address+output_file) as f:
            f[tree_name] = array_chunk

        print(f"Saved: {output_file}")
