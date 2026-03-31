import struct
import numpy as np


def build_W_1D(Lx, Ly, Nx, Ny, Nz, filename='z_coordinates.txt', dof=3, dtype=np.float64):
    """
    Build the 1D diagonal weight vector for the volume-weighted inner product.

    The domain is:
      - periodic and uniform in x and y
      - non-uniform in z, read from file

    The spatial flattening order matches the Fortran WRITE loop order:
      i (x) fastest, then j (y), then k (z)

    The final row ordering is assumed to be:
      [all U points, all V points, all W points]

    Returns
    -------
    W_1D : ndarray of shape (dof * Nx * Ny * Nz,)
    """

    if Nz < 2:
        raise ValueError("Nz must be at least 2.")

    try:
        z_coords = np.loadtxt(filename)
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find '{filename}'.")

    if z_coords.size != Nz:
        raise ValueError(
            f"Grid mismatch: expected Nz={Nz}, file has {z_coords.size} points."
        )

    dx = Lx / Nx
    dy = Ly / Ny

    total_nodes = Nx * Ny * Nz
    W_diag = np.zeros(total_nodes, dtype=dtype)

    idx = 0
    for k in range(Nz):
        if k == 0:
            dz = 0.5 * (z_coords[1] - z_coords[0])
        elif k == Nz - 1:
            dz = 0.5 * (z_coords[k] - z_coords[k - 1])
        else:
            dz = 0.5 * (z_coords[k + 1] - z_coords[k - 1])

        node_volume = dx * dy * dz

        for j in range(Ny):
            for i in range(Nx):
                W_diag[idx] = node_volume
                idx += 1

    domain_volume_calculated = np.sum(W_diag)
    domain_volume_actual = Lx * Ly * (z_coords[-1] - z_coords[0])

    print(f"Total weight sum : {domain_volume_calculated:.12e}")
    print(f"Actual volume    : {domain_volume_actual:.12e}")

    if np.isclose(domain_volume_calculated, domain_volume_actual, rtol=1e-10, atol=1e-12):
        print("SUCCESS: weight matrix matches physical domain volume.")
    else:
        print("WARNING: volumes do not match — check grid inputs.")

    W_1D = np.tile(W_diag, dof)
    return W_1D


def read_fortran_record(f, dtype):
    """
    Read one Fortran sequential unformatted record.
    Assumes 4-byte record markers before and after the data.
    """
    marker = f.read(4)
    if not marker:
        return None

    if len(marker) != 4:
        raise EOFError("Unexpected EOF while reading leading record marker.")

    nbytes = struct.unpack("<i", marker)[0]
    data = f.read(nbytes)

    if len(data) != nbytes:
        raise EOFError("Unexpected EOF while reading record payload.")

    end_marker_raw = f.read(4)
    if len(end_marker_raw) != 4:
        raise EOFError("Unexpected EOF while reading trailing record marker.")

    end_marker = struct.unpack("<i", end_marker_raw)[0]

    if end_marker != nbytes:
        raise ValueError(
            f"Record marker mismatch: start={nbytes}, end={end_marker}"
        )

    arr = np.frombuffer(data, dtype=dtype).copy()
    return arr


def build_Q_memmap(
    infile,
    nx,
    ny,
    nz,
    ntime,
    dtype=np.float32,
    records_per_timestep=3,
    out_memmap="Q_matrix.dat",
):
    """
    Reconstruct the full time-domain snapshot matrix Q from a Fortran binary file.

    Assumes each timestep is written as:
        record 1 = U field
        record 2 = V field
        record 3 = W field

    and each record contains nx*ny*nz values.

    Output
    ------
    Q : memmap of shape (3*nx*ny*nz, ntime)

    Row ordering:
        [all U points, all V points, all W points]
    Column ordering:
        one column per time step
    """

    ns_one = nx * ny * nz
    ns_total = records_per_timestep * ns_one

    Q = np.memmap(out_memmap, dtype=dtype, mode="w+", shape=(ns_total, ntime))

    with open(infile, "rb") as f:
        for t in range(ntime):
            blocks = []

            for r in range(records_per_timestep):
                rec = read_fortran_record(f, dtype)

                if rec is None:
                    raise EOFError(
                        f"Unexpected EOF at timestep {t + 1}, record {r + 1}"
                    )

                if rec.size != ns_one:
                    raise ValueError(
                        f"Record size mismatch at timestep {t + 1}, record {r + 1}: "
                        f"got {rec.size}, expected {ns_one}"
                    )

                blocks.append(rec)

            Q[:, t] = np.concatenate(blocks)

    Q.flush()
    return Q


def get_block_starts(ntime, block_size, overlap):
    """
    Return the starting indices of all complete blocks.
    """
    if block_size <= 1:
        raise ValueError("block_size must be at least 2.")
    if overlap < 0:
        raise ValueError("overlap must be non-negative.")
    if overlap >= block_size:
        raise ValueError("overlap must be smaller than block_size.")
    if block_size > ntime:
        raise ValueError("block_size cannot exceed ntime.")

    step = block_size - overlap
    starts = list(range(0, ntime - block_size + 1, step))
    return starts


def compute_spod_reduced_matrices(
    q_memmap_file,
    shape_q,
    W_1D,
    block_size,
    overlap,
    q_dtype=np.float32,
    real_dtype=np.float64,
    complex_dtype=np.complex128,
    subtract_block_mean=True,
    use_hann_window=True,
    out_file="R_freq.dat",
):
    """
    Compute the SPOD reduced matrices through step 8:

        1. split Q into overlapping time blocks
        2. optionally subtract block mean
        3. optionally apply Hann window
        4. FFT each block along time
        5. for each frequency, assemble X(f)
        6. compute R(f) = X(f)^* W X(f)

    Parameters
    ----------
    q_memmap_file : str
        Path to Q memmap file of shape (ns, ntime).

    shape_q : tuple
        Shape of Q as (ns, ntime).

    W_1D : ndarray
        1D diagonal weights of length ns.

    block_size : int
        Number of snapshots per block.

    overlap : int
        Number of overlapping snapshots between adjacent blocks.

    q_dtype : dtype
        Data type used to store Q.

    real_dtype : dtype
        Real working precision for preprocessing.

    complex_dtype : dtype
        Complex working precision for FFT and reduced matrices.

    subtract_block_mean : bool
        If True, subtract mean over time inside each block.

    use_hann_window : bool
        If True, apply Hann window in time before FFT.

    out_file : str
        Output memmap file for the frequency-wise reduced matrices.

    Returns
    -------
    R_freq : memmap
        Array of shape (nfreq, nblk, nblk), where
        R_freq[j] = X(f_j)^* W X(f_j)
    freq_values : ndarray
        Frequency bins corresponding to rFFT, in cycles per sample.
        Multiply by sampling rate if needed.
    starts : list
        Block starting indices.
    """

    ns, ntime = shape_q

    if W_1D.shape[0] != ns:
        raise ValueError(
            f"W_1D length mismatch: got {W_1D.shape[0]}, expected {ns}"
        )

    starts = get_block_starts(ntime, block_size, overlap)
    nblk = len(starts)
    nfreq = block_size // 2 + 1

    print(f"Number of blocks       : {nblk}")
    print(f"Block size             : {block_size}")
    print(f"Overlap                : {overlap}")
    print(f"Number of frequencies  : {nfreq}")

    Q = np.memmap(q_memmap_file, dtype=q_dtype, mode="r", shape=(ns, ntime))

    # Frequency-wise reduced matrices: one (nblk x nblk) matrix per frequency
    R_freq = np.memmap(
        out_file,
        dtype=complex_dtype,
        mode="w+",
        shape=(nfreq, nblk, nblk),
    )
    R_freq[:] = 0.0 + 0.0j

    if use_hann_window:
        window = np.hanning(block_size).astype(real_dtype)
    else:
        window = np.ones(block_size, dtype=real_dtype)

    # Build all block FFTs once, frequency by frequency
    # Xf_all[j] has shape (ns, nblk)
    Xf_all = [
        np.zeros((ns, nblk), dtype=complex_dtype)
        for _ in range(nfreq)
    ]

    for m, start in enumerate(starts):
        stop = start + block_size

        # Q_block shape: (ns, block_size)
        Q_block = np.asarray(Q[:, start:stop], dtype=real_dtype)

        if subtract_block_mean:
            Q_block = Q_block - Q_block.mean(axis=1, keepdims=True)

        Q_block = Q_block * window[np.newaxis, :]

        # FFT along time axis
        Qhat_block = np.fft.rfft(Q_block, axis=1).astype(complex_dtype)

        for j in range(nfreq):
            Xf_all[j][:, m] = Qhat_block[:, j]

        print(f"Processed FFT block {m + 1} / {nblk} (snapshots {start}:{stop})")

    # Step 8: compute R(f_j) = X(f_j)^* W X(f_j)
    for j in range(nfreq):
        Xf = Xf_all[j]                            # shape: (ns, nblk)
        WXf = W_1D[:, np.newaxis] * Xf           # shape: (ns, nblk)
        R_freq[j, :, :] = Xf.conj().T @ WXf      # shape: (nblk, nblk)

        print(f"Built weighted reduced matrix for frequency {j + 1} / {nfreq}")

    R_freq.flush()

    freq_values = np.fft.rfftfreq(block_size, d=1.0)

    return R_freq, freq_values, starts


if __name__ == "__main__":
    # =========================================================
    # USER INPUTS
    # =========================================================
    infile = "flattened_data.bin"

    nx = 192
    ny = 192
    nz = 129
    ntime = 400

    # Domain lengths
    # Replace Ly if needed
    Lx = 3.0 * np.pi
    Ly = 3.0 * np.pi

    # Binary data type
    # Use np.float32 if Fortran REAL was 4 bytes
    # Use np.float64 if the Fortran output used 8-byte reals
    q_dtype = np.float32

    # Working precisions
    w_dtype = np.float64
    real_work_dtype = np.float64
    complex_work_dtype = np.complex128

    # Raw binary layout
    records_per_timestep = 3
    dof = 3

    # SPOD block settings
    block_size = 50
    overlap = 25

    # Preprocessing flags
    subtract_block_mean = True
    use_hann_window = True

    # File names
    z_file = "z_coordinates.txt"
    q_memmap_file = "Q_matrix.dat"
    r_freq_file = "R_freq.dat"

    # =========================================================
    # DERIVED SIZES
    # =========================================================
    ns_one = nx * ny * nz
    ns_total = records_per_timestep * ns_one

    # =========================================================
    # STEP 1: BUILD W_1D
    # =========================================================
    print("\nBuilding W_1D ...")
    W_1D = build_W_1D(
        Lx=Lx,
        Ly=Ly,
        Nx=nx,
        Ny=ny,
        Nz=nz,
        filename=z_file,
        dof=dof,
        dtype=w_dtype,
    )
    print(f"W_1D shape: {W_1D.shape}")

    # =========================================================
    # STEP 2: BUILD FULL TIME-DOMAIN SNAPSHOT MATRIX Q
    # =========================================================
    print("\nBuilding full time-domain snapshot matrix Q ...")
    Q = build_Q_memmap(
        infile=infile,
        nx=nx,
        ny=ny,
        nz=nz,
        ntime=ntime,
        dtype=q_dtype,
        records_per_timestep=records_per_timestep,
        out_memmap=q_memmap_file,
    )
    print(f"Q stored in: {q_memmap_file}")
    print(f"Q shape: ({ns_total}, {ntime})")

    # =========================================================
    # STEPS 3 THROUGH 8: BLOCKING, FFT, ASSEMBLY, WEIGHTED MATRICES
    # =========================================================
    print("\nComputing SPOD reduced matrices through step 8 ...")
    R_freq, freq_values, starts = compute_spod_reduced_matrices(
        q_memmap_file=q_memmap_file,
        shape_q=(ns_total, ntime),
        W_1D=W_1D,
        block_size=block_size,
        overlap=overlap,
        q_dtype=q_dtype,
        real_dtype=real_work_dtype,
        complex_dtype=complex_work_dtype,
        subtract_block_mean=subtract_block_mean,
        use_hann_window=use_hann_window,
        out_file=r_freq_file,
    )

    print("\nDone through step 8.")
    print(f"Q memmap file          : {q_memmap_file}")
    print(f"R_freq memmap file     : {r_freq_file}")
    print(f"Block starts           : {starts}")
    print(f"Number of frequencies  : {len(freq_values)}")
    print(f"R_freq shape           : {R_freq.shape}")