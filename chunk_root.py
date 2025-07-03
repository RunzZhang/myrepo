import uproot
import awkward as ak
import numpy as np
base_address = "/data/runzezhang/result/TN_sims_D/"
root_name = "dmx_Cf_1E7.root"
input_file = base_address+root_name

tree_name = "tree"  # Replace with your tree name
columns = None  # or list of specific branches if you want to filter
num_parts = 20
def to_fixed_str(array, length=32):
    return np.array(array, dtype=f"S{length}")  # e.g. 32-byte strings

def bytes_to_unicode(array, length=32):
    # Decode from bytes to Unicode string (U)
    return np.char.decode(array, encoding='utf-8').astype(f"<U{length}")

def to_fixed_unicode(array, length=32):
    # Convert to NumPy array if not already
    np_array = np.asarray(array)

    # If it's object type or bytes, decode manually
    # if np_array.dtype.kind in {'O', 'S'}:
    #     decoded = np.char.decode(np_array, encoding='utf-8', errors='replace')
    #     return decoded.astype(f"<U{length}")
    # elif np_array.dtype.kind == 'U':
    #     return np_array.astype(f"<U{length}")
    # else:
    #     raise TypeError(f"Unexpected dtype for string field: {np_array.dtype}")

    # return np.array(array, dtype=f"<U{length}")
    return np.array(ak.to_list(array), dtype=f"<U{length}")


def clean_array_chunk(array_chunk, str_length=64):
    cleaned = {}
    for field in array_chunk.fields:
        data = array_chunk[field]

        # Convert to plain list (safe)
        values = ak.to_list(data)

        # Decide dtype
        if isinstance(values[0], str):
            cleaned[field] = np.array(values, dtype=f"<U{str_length}")
        elif isinstance(values[0], (bytes, np.bytes_)):
            decoded = [v.decode('utf-8', errors='replace') for v in values]
            cleaned[field] = np.array(decoded, dtype=f"<U{str_length}")
        else:
            cleaned[field] = np.array(values)

    return cleaned  # This is a dict of NumPy arrays

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
        # array_chunk = ak.with_field(array_chunk, to_fixed_unicode(array_chunk["Process"]), "Process")
        # array_chunk = ak.with_field(array_chunk, to_fixed_unicode(array_chunk["name"]), "name")
        # array_chunk = ak.with_field(array_chunk, to_fixed_unicode(array_chunk["Volume"]), "Volume")
        array_cleaned = clean_array_chunk(array_chunk)


        # Save to new root file
        output_file = f"dmx_Cf_1e7_part{i}.root"
        with uproot.recreate(base_address+output_file) as f:
            f[tree_name] = array_chunk

        print(f"Saved: {output_file}")
