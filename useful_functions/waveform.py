import numpy as np
from useful_functions import read_and_extract_parameters

def read_vectors(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        vectors = []
        header = lines[:4]
        for e, line in enumerate(lines):
            if e > 3:
                try:
                    vectors.append([float(x) for x in line.strip().split(',')])
                except ValueError:
                    continue
        return header, np.array(vectors)

def find_first_non_zero_vector(vectors):
    for vec in vectors:
        if not np.allclose(vec, 0):
            return vec
    return None

def rotation_matrix_from_vectors(vec1, vec2):
    """ 
    Find the rotation matrix that aligns vec1 to vec2 
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def transform_vectors(vectors, rotation_matrix):
    return np.dot(vectors, rotation_matrix.T)

def write_vectors(header, file_path, vectors):
    with open(file_path, 'w') as file:
        for h in header:
            file.write(h)
        for vec in vectors:
            file.write(' '.join(map(str, vec)) + '\n')

def write_params(file_path, bs, vecs):
    with open(file_path, 'w') as file:
        for b in bs :
            for vec in vecs:
                # write vec and b 
                file.write(f'{vec[0]} {vec[1]} {vec[2]} {b}\n')


def create_waveform(input_file,  target_vector, scale):
    header, vectors = read_vectors(input_file)
    first_non_zero_vector = find_first_non_zero_vector(vectors)
    if first_non_zero_vector is None:
        print("No non-zero vector found in the file.")
        return
    
    rotation_matrix = rotation_matrix_from_vectors(first_non_zero_vector, target_vector)
    transformed_vectors = transform_vectors(vectors, rotation_matrix)
    #rescale all vectors
    transformed_vectors *= scale

    return header, transformed_vectors
    

if __name__ == "__main__":
    # Example usage
    input_file = '/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/full_waveform_flipped.txt'
    bvalues = [200, 1000]
    scheme_file = '/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/PGSE_21_dir_12_b.scheme'

    _,_, _, _, vectors = read_and_extract_parameters(scheme_file)

    # Convert the 2D array to a structured array for easy comparison
    dtype = np.dtype((np.void, vectors.dtype.itemsize * vectors.shape[1]))
    structured_arr = np.ascontiguousarray(vectors).view(dtype)
    
    # Use numpy.unique to find unique rows
    _, idx = np.unique(structured_arr, return_index=True)
    
    # Return the unique rows
    vectors = vectors[idx]
    
    all_waveforms = []
    for b in bvalues:
        if b == 200:
            scale = 0.022447536757108
        elif b == 1000:
            scale = 0.050194218116318
        else:
            #error
            assert False
        for e,vector in enumerate(vectors):
         
            header, waveforms = create_waveform(input_file, vector, scale)
            all_waveforms.extend(waveforms)
    
    output_file = f'/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/iso_waveform.txt'
    write_vectors(header, output_file, all_waveforms)
    params_path = '/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/iso_waveform_vec_b.txt'
    write_params(params_path, bvalues, vectors)
    print(f"Transformation applied and saved to {output_file}")