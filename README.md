\documentclass[10pt,a4paper]{article}
\usepackage[margin=0.75in]{geometry}
\usepackage{amsmath, amssymb, bm}
\usepackage{titlesec}
\usepackage{enumitem}
\usepackage{multicol}
\usepackage{tikz}
\usepackage{booktabs}
\usepackage{longtable}
\usetikzlibrary{shapes.geometric, arrows, positioning}

\titleformat{\section}{\large\bfseries}{\thesection}{0.6em}{}[\titlerule]
\titlespacing*{\section}{0pt}{1.5ex plus 1ex minus .2ex}{0.8ex plus .2ex}

\begin{document}

\begin{center}
    {\LARGE \textbf{Spectral Proper Orthogonal Decomposition\\
                    of Wind-over-Waves DNS Data}} \\
    \vspace{2mm}
    \textbf{Sudin Bhuju Shrestha} \\
    \textit{April 2026}
\end{center}

\vspace{2mm}

%==========================================================
\section{Introduction}
\label{sec:intro}
%==========================================================

\subsection{Motivation}
\label{sec:motivation}

Standard Proper Orthogonal Decomposition (POD) identifies the most 
energetic spatial structures of a flow field but ignores temporal 
evolution; modes that mix multiple physical mechanisms operating at 
different frequencies are not separated by the decomposition. 
Spectral POD (SPOD) addresses this limitation by retaining both 
spatial and temporal information. Applying a temporal Fourier 
transform to the space--time correlation tensor yields the 
Cross-Spectral Density (CSD) tensor $\hat{C}(\bm{x}, \bm{x}', f)$, 
whose eigendecomposition produces modes that are coherent in both 
space and time. SPOD therefore isolates physically meaningful 
phenomena --- such as shedding vortices, traveling waves, or 
wave-coherent boundary-layer structures --- at specific frequencies.

This report applies SPOD to a Direct Numerical Simulation (DNS) 
dataset of a wind-over-waves boundary layer. The grid contains 
$N_x \times N_y \times N_z = 192 \times 192 \times 65$ points, and 
the snapshot ensemble comprises $500$ temporal samples covering the 
streamwise, spanwise, and wall-normal velocity components 
$(U, V, W)$. The objective is to identify the dominant coherent 
structures of the boundary layer, quantify their wall-normal 
distribution, and characterise the transition from coherent to 
incoherent dynamics across the resolved frequency band.

\subsection{Welch's Method for Cross-Spectral Density Estimation}
\label{sec:welch}

A reliable estimate of the CSD requires an ensemble of statistically 
independent realisations. Since the DNS produces a single continuous 
time series, Welch's method is employed to construct this ensemble. 
The full record is partitioned into $N_{\mathrm{blk}}$ overlapping 
blocks, each of which is independently Fourier-transformed in time. 
This procedure provides a controlled compromise between frequency 
resolution and statistical noise reduction: longer blocks improve 
resolution at the cost of fewer realisations, while more numerous 
short blocks reduce variance at the cost of resolution. For a 
selected frequency $f_j$, the Fourier coefficients across all blocks 
are gathered into the frequency-domain snapshot matrix
\begin{equation*}
    X(f_j) = 
    \left[\hat{\bm{q}}^{(1)}(f_j),\,
          \hat{\bm{q}}^{(2)}(f_j),\,
          \dots,\,
          \hat{\bm{q}}^{(N_{\mathrm{blk}})}(f_j)\right] 
    \in \mathbb{C}^{3N \times N_{\mathrm{blk}}}.
\end{equation*}

\subsection{The Spatial Weight Matrix}
\label{sec:weight_motivation}

The weight matrix $W$ defines the inner product under which the 
spatial modes are optimal and orthogonal. In an incompressible-flow 
SPOD analysis, $W$ acts as a discrete spatial quadrature operator 
that accounts for the volume of influence of each grid point. This 
correction is essential whenever the computational mesh is 
non-uniform: in the present DNS, the wall-normal direction is 
clustered near the surface to resolve the boundary-layer shear, so 
that without spatial weighting the densely-packed near-surface nodes 
would dominate the kinetic-energy calculation purely due to grid 
density rather than physical content. The construction of $W$ for 
the present grid is detailed in Section~\ref{sec:weight_compute}.

%==========================================================
\section{Theoretical Formulation}
\label{sec:theory}
%==========================================================

This section presents the weighted Method of Snapshots, which 
forms the mathematical backbone of the SPOD pipeline. The 
derivation tracks how the grid weight matrix $W$ enters both the 
reduced eigenvalue problem and the final spatial mode 
reconstruction.

\subsection{Weighted Method of Snapshots}
\label{sec:method_of_snapshots}

The theoretical weighted SPOD eigenvalue problem at a specific 
frequency $f_j$ is
\begin{equation}
    \hat{C} W \hat{\Phi} = \hat{\Phi} \hat{\Lambda}.
    \label{eq:spod_eigenvalue}
\end{equation}
Substituting the empirical CSD estimate $\hat{C} = \frac{1}{N_{\mathrm{blk}}} X X^*$ yields
\begin{equation}
    \frac{1}{N_{\mathrm{blk}}} X X^* W \hat{\Phi} = \hat{\Phi} \hat{\Lambda}.
    \label{eq:spod_substituted}
\end{equation}
Because $X X^* W$ is an asymmetric matrix of size $3N \times 3N$, 
this system cannot be solved efficiently by direct decomposition. 
Multiplying both sides on the left by $X^* W$ gives
\begin{equation}
    X^* W \left( \frac{1}{N_{\mathrm{blk}}} X X^* W \hat{\Phi} \right) 
    = X^* W \hat{\Phi} \hat{\Lambda},
\end{equation}
and using associativity,
\begin{equation}
    \frac{1}{N_{\mathrm{blk}}} \left( X^* W X \right) 
    \left( X^* W \hat{\Phi} \right) 
    = \left( X^* W \hat{\Phi} \right) \hat{\Lambda}.
\end{equation}
Defining the reduced correlation matrix 
$R \equiv \frac{1}{N_{\mathrm{blk}}} X^* W X$ and the reduced 
eigenvectors $\Psi \equiv X^* W \hat{\Phi}$, this becomes
\begin{equation}
    R \Psi = \Psi \hat{\Lambda}.
    \label{eq:reduced_eigenvalue}
\end{equation}

The reduced system is computationally tractable because $R$ is 
only $N_{\mathrm{blk}} \times N_{\mathrm{blk}}$. The physical 
influence of the grid weights $W$ is embedded inside $R$. Solving 
Equation~\eqref{eq:reduced_eigenvalue} yields the modal energies 
$\hat{\Lambda}$ and the orthonormal reduced eigenvectors $\Psi$ 
(satisfying $\Psi^* \Psi = I$).

\subsubsection*{Spatial mode reconstruction}

Because $X^* W$ is rectangular, the definition $\Psi = X^* W \hat{\Phi}$ 
cannot be inverted directly to recover the spatial modes. 
Following Sirovich's snapshot-method assumption, the spatial modes 
are constructed as a linear combination of the original frequency 
snapshots, weighted by the reduced eigenvectors:
\begin{equation}
    \hat{\Phi}_{\mathrm{raw}} = X \Psi.
    \label{eq:phi_raw}
\end{equation}

\subsubsection*{Orthonormality scaling}

The final spatial modes must satisfy the weighted discrete 
orthonormality condition $\hat{\Phi}^* W \hat{\Phi} = I$. Testing 
the raw modes,
\begin{equation}
    \hat{\Phi}_{\mathrm{raw}}^* W \hat{\Phi}_{\mathrm{raw}} 
    = (X \Psi)^* W (X \Psi) 
    = \Psi^* (X^* W X) \Psi.
\end{equation}
From Equation~\eqref{eq:reduced_eigenvalue}, 
$X^* W X \, \Psi = N_{\mathrm{blk}} \Psi \hat{\Lambda}$, giving
\begin{equation}
    \hat{\Phi}_{\mathrm{raw}}^* W \hat{\Phi}_{\mathrm{raw}} 
    = N_{\mathrm{blk}} \Psi^* \Psi \hat{\Lambda} 
    = N_{\mathrm{blk}} \hat{\Lambda}.
\end{equation}
The squared length of the raw modes is therefore scaled by 
$N_{\mathrm{blk}} \hat{\Lambda}$. Because the orthonormality 
condition is quadratic in $\hat{\Phi}$, the correctly scaled SPOD 
modes are
\begin{equation}
    \hat{\Phi} 
    = \frac{1}{\sqrt{N_{\mathrm{blk}} \hat{\Lambda}}} \, X \Psi.
    \label{eq:phi_normalized}
\end{equation}
Equation~\eqref{eq:reduced_eigenvalue} delivers the modal energies 
$\hat{\Lambda}$, and Equation~\eqref{eq:phi_normalized} delivers 
the spatial modes $\hat{\Phi}$.

\subsection{Procedure Summary}
\label{sec:procedure_summary}

The complete SPOD procedure derived above is summarised 
schematically in Figure~\ref{fig:flowchart}.

\begin{figure}[hbt!]
\centering
\begin{tikzpicture}[node distance=1.25cm, auto]

\tikzstyle{data} = [trapezium, trapezium left angle=70, trapezium right angle=110, minimum width=4cm, minimum height=1cm, text centered, draw=black, fill=blue!10]
\tikzstyle{process} = [rectangle, rounded corners, minimum width=6cm, minimum height=1cm, text centered, draw=black, fill=gray!10]
\tikzstyle{math} = [rectangle, minimum width=6cm, minimum height=1cm, text centered, draw=black, fill=green!10]
\tikzstyle{arrow} = [thick,->,>=stealth]

\node (raw) [data] {Raw Spatial--Temporal Data: $q(\mathbf{x}, t)$};
\node (fluct) [process, below of=raw] {Compute Fluctuations: $q'(\mathbf{x}, t) = q - \bar{q}$};
\node (blocks) [process, below of=fluct] {Segment into $N_{\mathrm{blk}}$ overlapping blocks};
\node (window) [process, below of=blocks] {Apply temporal window function (Hann)};
\node (fft) [process, below of=window] {Perform temporal FFT on each block};
\node (matrix) [process, below of=fft] {Assemble Frequency Data Matrix: $X(f_j)$};
\node (corr) [math, below of=matrix] {Form Correlation Matrix: $R(f_j) = \frac{1}{N_{\mathrm{blk}}} X^*(f_j) W X(f_j)$};
\node (eigen) [math, below of=corr] {Solve Eigenproblem: $R(f_j) \Psi(f_j) = \Psi(f_j) \hat{\Lambda}(f_j)$};
\node (reconstruct) [math, below of=eigen]
  {Reconstruct Modes: $\hat{\Phi}(f_j) = \left(N_{\mathrm{blk}} \hat{\Lambda}(f_j)\right)^{-1/2} X(f_j)\, \Psi(f_j)$};
\node (output) [data, below of=reconstruct] {Output: Modes $\hat{\Phi}(f_j)$ \& Energies $\hat{\Lambda}(f_j)$};

\draw [arrow] (raw) -- (fluct);
\draw [arrow] (fluct) -- (blocks);
\draw [arrow] (blocks) -- (window);
\draw [arrow] (window) -- (fft);
\draw [arrow] (fft) -- (matrix);
\draw [arrow] (matrix) -- (corr);
\draw [arrow] (corr) -- (eigen);
\draw [arrow] (eigen) -- (reconstruct);
\draw [arrow] (reconstruct) -- (output);

\end{tikzpicture}
\caption{SPOD procedure flowchart. Blue trapezoids denote 
data inputs and outputs; grey blocks denote preprocessing steps; 
green blocks denote the core mathematical operations of the weighted 
Method of Snapshots.}
\label{fig:flowchart}
\end{figure}

%==========================================================
\section{Computational Implementation}
\label{sec:implementation}
%==========================================================

This section bridges the theoretical formulation of 
Section~\ref{sec:theory} with a memory-efficient discrete 
implementation. The pipeline is structured around four phases: 
data ingestion and matrix setup, preprocessing (windowing and 
fluctuation extraction), spectral estimation, and mode 
reconstruction.

\subsection{Data Ingestion and Matrix Setup}
\label{sec:data_ingestion}

The raw DNS data, originally generated by Fortran using unformatted 
stream access (\texttt{ACCESS='STREAM'}), exists as a continuous 
sequence of binary bytes without record markers. Because the SPOD 
algorithm requires rapid, non-sequential access to overlapping time 
segments, the data must be structurally rearranged. The algorithm 
systematically reads the raw stream block by block, identifying 
exactly $N_x \times N_y \times N_z$ scalar values for each of the 
three velocity components per timestep. These sequential blocks are 
concatenated to build the structured 2D time-domain snapshot matrix 
$Q$, where rows represent flattened spatial coordinates and columns 
represent discrete time snapshots. To manage the data footprint and 
allow rapid access to time blocks, the rearranged matrix $Q$ is 
written directly to disk as a memory-mapped array via the 
\texttt{build\_Q\_memmap} routine, circumventing system memory 
limitations during processing.

\subsection{Spatial Weight Matrix Computation}
\label{sec:weight_compute}

The discrete weight matrix $W$ defines the inner product under 
which the SPOD eigenvalue problem~\eqref{eq:spod_eigenvalue} is 
posed. Physically, $W$ acts as a numerical quadrature rule that 
accounts for the varying spatial measure of each grid point.

For the present computational domain, the horizontal grid spacing 
is uniform ($\Delta x = L_x / N_x$ and $\Delta y = L_y / N_y$), 
while the vertical coordinates $z_k$ are strictly increasing but 
non-uniform. The discrete vertical weight at index $k$ is computed 
using the trapezoidal rule:
\begin{equation}
    \Delta z_k = \tfrac{1}{2} (z_{k+1} - z_{k-1}),
    \label{eq:dz_trap}
\end{equation}
with one-sided differences applied at the domain boundaries. The 
local spatial weight for any node is the three-dimensional control 
volume $\Delta V_{ijk} = \Delta x \, \Delta y \, \Delta z_k$.

Algorithmic consistency is verified by ensuring that the sum of all 
nodal weights exactly reproduces the physical domain volume:
\begin{equation}
    \sum_{k=0}^{N_z-1} \sum_{j=0}^{N_y-1} \sum_{i=0}^{N_x-1} 
    \Delta x \, \Delta y \, \Delta z_k 
    = L_x L_y \, (z_{\max} - z_{\min}).
    \label{eq:weight_check}
\end{equation}
For the present grid, the computed weight sum matches the analytical 
value $9\pi^2 \approx 88.83$ to machine precision.

Because the state vector consists of three dimensionally identical 
velocity components, no specialised variable-scaling is required. 
The full weight matrix $W$ is constructed as a 1D diagonal array by 
tiling the localised grid volumes across the three velocity 
components, ensuring that every component is scaled identically by 
the physical mesh geometry.

\subsection{Block Segmentation and Fluctuation Extraction}
\label{sec:block_segmentation}

To estimate a robust CSD, the snapshot matrix $Q$ is partitioned 
into an ensemble of shorter temporal blocks. If the blocks were 
extracted sequentially without overlap (as in the Bartlett method), 
the abrupt truncation at the block boundaries would introduce 
artificial high-frequency noise into the Fourier transform --- 
spectral leakage. To mitigate this, a Hann window is applied to 
each block, smoothly tapering the velocity signal to zero at the 
beginning and end of the segment and enforcing strict periodicity.

The Hann window introduces a secondary complication: by forcing the 
boundaries to zero, it discards physical content near the edges of 
every block. To compensate, the algorithm uses Welch's method with 
50\% block overlap, ensuring that turbulence attenuated at the edge 
of one block sits at the unattenuated peak of an adjacent one. The 
overlap also increases the number of available block realisations, 
raising the statistical degrees of freedom and reducing the variance 
of the spectral estimate. Because true turbulent structures are 
coherent across time, ensemble-averaging across many overlapping 
blocks suppresses random noise while preserving distinct energy 
peaks of the underlying physical phenomena.

The matrix is loaded into memory piece by piece by selectively 
extracting block-specific columns from the memmap. Once a block is 
loaded, the temporal mean of that specific block is calculated and 
subtracted from every row, applying a Reynolds decomposition to 
convert raw velocity data into purely turbulent fluctuation data:
\begin{equation}
    q'(\mathbf{x}, t) 
    = q(\mathbf{x}, t) - \bar{q}_{\mathrm{block}}(\mathbf{x}).
    \label{eq:reynolds}
\end{equation}

Performing the mean subtraction at this stage --- rather than 
having the Fortran solver export pre-computed fluctuations --- is a 
deliberate architectural choice. Configuring the DNS to export the 
full instantaneous velocity field preserves the complete physical 
state of the flow and maximises the dataset's flexibility, allowing 
the same binary output to be reused for different analyses without 
re-running the simulation. The SPOD pipeline therefore takes 
responsibility for isolating the fluctuations.

From a fluid-dynamics perspective, SPOD is designed to extract 
dynamic, coherent turbulent structures interacting at specific 
non-zero frequencies. If the stationary mean flow $\bar{q}$ were 
retained, its kinetic energy would dominate the decomposition and 
manifest as an overwhelming zero-frequency mode. By subtracting the 
block mean, the algorithm isolates fluctuations relative to the 
instantaneous baseline velocity profile and directs the 
decomposition toward transient eddy structures at non-zero 
frequencies. The residual energy visible at $f = 0$ reflects 
block-to-block variations in the mean flow that persist despite 
this preprocessing --- a feature expected when the baseline state 
is not strictly stationary.

\subsection{Spectral Transformation and Frequency Matrix Assembly}
\label{sec:fft_assembly}

Following preprocessing, the algorithm transforms the data into the 
frequency domain. For every overlapping block, a discrete real FFT 
is applied along the temporal dimension. To preserve physical 
energy, the resulting Fourier coefficients are scaled by the factor
\begin{equation}
    \kappa = \sqrt{\dfrac{\Delta t}{\sum_{n=1}^{N_{\mathrm{blk}}} w_n^2}},
    \label{eq:psd_scaling}
\end{equation}
where $w_n$ is the amplitude of the Hann window at sample $n$ and 
$\Delta t$ is the snapshot time spacing. This scaling compensates 
for the variance attenuation introduced by the windowing and 
ensures that the final eigenvalues $\hat{\Lambda}(f)$ represent 
power spectral density with correct physical dimensions.

In classical SPOD theory, applying the FFT to all blocks yields a 
three-dimensional tensor of size 
$N_{\mathrm{space}} \times N_{\mathrm{blk}} \times N_f$. Holding 
this tensor in memory is computationally prohibitive for 
high-fidelity DNS datasets. The present implementation evaluates 
the decomposition one frequency at a time: for a targeted discrete 
frequency $f_j$, the algorithm extracts the corresponding complex 
Fourier coefficients from every block and concatenates them into 
the localised 2D frequency-domain snapshot matrix 
$X(f_j) \in \mathbb{C}^{N_{\mathrm{space}} \times N_{\mathrm{blk}}}$. 
This isolated matrix contains all the spatial and statistical 
information required for the reduced eigenvalue problem at the 
chosen frequency.

\subsection{Assembly of the Reduced Correlation Matrix}
\label{sec:R_assembly}

The pipeline is structured to solve the reduced eigenvalue problem 
of Equation~\eqref{eq:reduced_eigenvalue}. While the theoretical 
construction of $R$ assumes a monolithic matrix product 
$X^* W X$, computing this simultaneously over the full spatial 
domain exceeds standard memory limits.

To circumvent this, the algorithm exploits the distributive 
property of matrix multiplication. The spatial domain is divided 
into manageable segments (chunks of $20{,}000$ grid points). For 
each spatial chunk, the algorithm extracts the corresponding rows 
of $X(f_j)$, multiplies them by their localised spatial weights 
from $W$, and computes the partial inner product. These partial 
matrices are cumulatively summed across all spatial chunks and 
finally ensemble-averaged by the block count $N_{\mathrm{blk}}$. 
This chunking strategy yields the exact $R$ matrix without ever 
loading the full $X(f_j)$ into memory simultaneously.

While $R$ is theoretically guaranteed to be Hermitian, the chunked 
floating-point operations introduce microscopic numerical roundoff. 
To ensure strict numerical stability for the eigenvalue solver, 
Hermitian symmetry is explicitly enforced on the final matrix by 
averaging it with its conjugate transpose, $R \leftarrow \tfrac{1}{2}(R + R^*)$, 
before decomposition.

\subsection{Reduced Eigenvalue Decomposition}
\label{sec:eigendecomposition}

Once $R(f_j)$ is constructed and stabilised, the algorithm solves 
the Hermitian eigenvalue problem of 
Equation~\eqref{eq:reduced_eigenvalue}. An optimised linear-algebra 
solver extracts the discrete modal energies $\hat{\Lambda}$ and the 
corresponding orthonormal reduced eigenvectors $\Psi$. Because 
SPOD's primary objective is to isolate the most dominant coherent 
flow structures, the resulting eigenvalues (and associated 
eigenvectors) are sorted in strictly descending order, so that the 
first mode $\hat{\Phi}_1$ is mathematically guaranteed to represent 
the most energetic phenomenon at the targeted frequency $f_j$.

A numerical-tolerance filter is applied to the eigenvalues. While 
true physical kinetic energy must be strictly non-negative, the 
floating-point operations involved in building $R$ can produce 
microscopic non-physical negative eigenvalues at the level of 
machine epsilon ($\sim 10^{-16}$). The algorithm identifies these 
numerical artifacts and forces them to zero, ensuring physical 
validity before the spatial reconstruction phase.

\subsection{Spatial Mode Reconstruction and Verification}
\label{sec:mode_reconstruction}

With $\hat{\Lambda}$ and $\Psi$ extracted, the pipeline reconstructs 
the 3D spatial modes $\hat{\Phi}$ via 
Equation~\eqref{eq:phi_normalized}. Just as for the correlation 
matrix assembly, computing 
$\hat{\Phi} = (N_{\mathrm{blk}} \hat{\Lambda})^{-1/2} \, X \Psi$ 
directly exceeds available memory. The pipeline therefore reads the 
frequency matrix $X(f_j)$ in spatial segments, projects each 
segment onto the scaled eigenvectors, and streams the resulting 
spatial modes directly to disk.

As a final validation step, the orthonormality of the reconstructed 
modes is verified by a localised, chunked computation of 
$\hat{\Phi}^* W \hat{\Phi}$ and comparison to the identity matrix. 
The maximum residual error of $\sim 10^{-15}$ confirms that the 
extracted modes are physically independent, orthogonal, and 
correctly scaled with respect to the non-uniform grid volume.

\subsection{Eigenvalues as Power Spectral Density}
\label{sec:psd_interpretation}

The eigenvalues $\hat{\Lambda}_k(f_j)$ computed at each frequency 
represent the contribution of the $k$-th spatial mode to the 
kinetic-energy spectrum at that frequency. Because the Fourier 
coefficients were scaled by $\kappa$ 
(Equation~\eqref{eq:psd_scaling}), the eigenvalues inherit the 
physical dimension of power spectral density: energy per unit 
frequency.

The sum $\sum_{k=1}^{N_{\mathrm{blk}}} \hat{\Lambda}_k(f_j)$ yields 
the total kinetic energy density at frequency $f_j$. Integrating 
across all frequencies recovers the total turbulent kinetic energy 
of the flow:
\begin{equation}
    \mathrm{TKE} 
    = \sum_{f_j} \sum_{k=1}^{N_{\mathrm{blk}}} 
      \hat{\Lambda}_k(f_j) \, \Delta f,
    \label{eq:tke}
\end{equation}
where $\Delta f$ is the frequency resolution. The SPOD spectrum 
therefore quantifies the distribution of kinetic energy across both 
spatial scales (mode index $k$) and temporal scales (frequency $f$), 
and the dominant eigenvalue at the peak frequency identifies both 
the location of the strongest coherent structure in the spectrum 
and its energetic importance relative to the turbulent cascade.

\subsection{Full Spectrum Execution}
\label{sec:full_spectrum_execution}

With the single-frequency pipeline validated, the algorithm 
iterates over the full set of discrete frequencies returned by the 
real FFT. The physical frequency axis is determined by the snapshot 
sampling rate $\Delta t$.

To ensure stable execution without exceeding memory limits or 
encountering operating-system file-lock conflicts, the pipeline 
explicitly purges memory-mapped file handles and triggers Python 
garbage collection at the start of each frequency iteration. The 
spatial orthonormality verification 
($\hat{\Phi}^* W \hat{\Phi} = I$) is restricted to the initial 
frequency step to minimise computational overhead, since the 
projection's mathematical structure is invariant across 
frequencies. As the pipeline iterates through the spectrum, the 
localised 3D spatial modes $\hat{\Phi}$ are written directly to 
disk while the modal energies $\hat{\Lambda}$ are aggregated into a 
global data structure. This aggregation yields the complete SPOD 
energy spectrum, which serves as the primary diagnostic tool for 
identifying coherent turbulent structures and their characteristic 
frequencies.

%==========================================================
\section{Results}
\label{sec:results}
%==========================================================

This section presents the SPOD analysis of the wind-over-waves DNS 
dataset. The discussion proceeds from a global characterisation of 
the eigenvalue spectrum through quantitative diagnostics of mode 
dominance, to the spatial structure of the leading modes and their 
component-energy distribution. The pipeline produced 26 frequency 
bins from $N_{\mathrm{blk}} = 19$ overlapping blocks, with frequency 
resolution $\Delta f = 0.02$ and Nyquist frequency $f_{\mathrm{Ny}} = 0.5$.

\subsection{Spectral Energy Distribution and Mode Dominance}
\label{sec:spectral_summary}

The modal energy distribution was analysed across the discrete 
frequency bins using three diagnostic ratios: the leading-mode 
energy $\hat{\Lambda}_1$, the mode-1 coherence ratio 
$\hat{\Lambda}_1 / \sum_k \hat{\Lambda}_k$, and the rank separation 
$\hat{\Lambda}_1 / \hat{\Lambda}_2$. The total summed modal energy 
across the spectrum is $\sum_f \sum_k \hat{\Lambda}_k(f) = 2.02 \times 10^{1}$. 
As shown in Table~\ref{tab:spod_top_freqs}, the spectrum exhibits a 
sharp global peak at $f = 0.02$.

\begin{table}[hbt!]
\centering
\caption{Top five frequencies ranked by leading mode energy 
$\hat{\Lambda}_1$. The $f = 0$ entry reflects residual block-to-block 
variation in the mean velocity field rather than a true 
zero-frequency mode, since each block undergoes independent mean 
subtraction prior to the FFT.}
\label{tab:spod_top_freqs}
\begin{tabular}{@{}ccccccc@{}}
\toprule
\textbf{Rank} & \textbf{$f$} & $\sum \hat{\Lambda}_k$ &
$\hat{\Lambda}_1$ &
$\hat{\Lambda}_1 / \sum \hat{\Lambda}_k$ &
$\hat{\Lambda}_1 / \hat{\Lambda}_2$ &
\textbf{Modes $\geq 90\%$} \\
\midrule
1 & 0.0200 & 3.474 & 1.178 & 0.339 & 1.97 & 10 \\
2 & 0.0400 & 3.079 & 0.466 & 0.151 & 1.14 & 14 \\
3 & 0.0000 & 1.655 & 0.359 & 0.217 & 1.21 & 8  \\
4 & 0.0600 & 2.256 & 0.261 & 0.116 & 1.19 & 16 \\
5 & 0.0800 & 1.719 & 0.254 & 0.148 & 1.53 & 16 \\
\bottomrule
\end{tabular}
\end{table}

At the peak frequency $f = 0.02$, the leading mode captures 
$33.9\%$ of the local spectral energy, and the rank separation 
$\hat{\Lambda}_1 / \hat{\Lambda}_2 = 1.97$ is just below the 
$\geq 2.0$ threshold conventionally used to define rank-1 dominance. 
The strongest separation in the entire spectrum occurs at $f = 0.12$, 
with $\hat{\Lambda}_1 / \hat{\Lambda}_2 = 2.09$, suggesting that 
this frequency hosts a particularly distinct coherent structure 
despite its lower absolute energy. These two frequencies, together 
with $f = 0.06$ (the broader coherent regime) and $f = 0.40$ (a 
turbulent reference), define the four target frequencies analysed 
in detail in subsequent sections.

The mode-1 coherence ratio identifies a clear coherent--incoherent 
transition. With $N_{\mathrm{blk}} = 19$ blocks, the incoherent 
baseline is $1/N_{\mathrm{blk}} = 0.053$ and the conventional 
coherence threshold is $2/N_{\mathrm{blk}} = 0.105$. The ratio 
$\hat{\Lambda}_1 / \sum_k \hat{\Lambda}_k$ exceeds this threshold 
only for $f \leq 0.12$, with the first crossing below it occurring 
at $f = 0.14$. Beyond this frequency, the energy is distributed 
nearly uniformly across all blocks, and the number of modes 
required to capture $90\%$ of the variance saturates at 16--17 --- 
close to the maximum of $N_{\mathrm{blk}} = 19$. This is the 
spectral signature of broadband, spatially uncorrelated turbulent 
fluctuations and defines the upper edge of the coherent regime at 
$f \approx 0.14$.

\subsection{SPOD Eigenvalue Spectrum}
\label{sec:spod_spectrum}

Figure~\ref{fig:spod_spectrum} presents the complete SPOD eigenvalue 
spectrum across all 26 resolved frequencies, displayed on both 
linear-frequency/log-energy and log--log axes.

\begin{figure}[hbt!]
    \centering
    \includegraphics[width=\textwidth]{figures/spod_spectrum.png}
    \caption{SPOD eigenvalue spectrum across all resolved frequencies. 
    Each curve corresponds to a single mode rank (modes 1--19, 
    darkest to lightest); the dashed black line is the total energy 
    per frequency $\sum_k \hat{\Lambda}_k$. Left: log-energy versus 
    linear frequency, with the global peak at $f = 0.02$ (red dotted 
    line). Right: log--log axes, with a dash-dotted $f^{-5/3}$ 
    reference line.}
    \label{fig:spod_spectrum}
\end{figure}

The spectrum reveals three distinct frequency regimes. At low 
frequencies ($f \leq 0.04$), the leading mode is well separated 
from the remaining modes, with the global energy peak at $f = 0.02$ 
where $\hat{\Lambda}_1 = 1.178$ and the coherence ratio reaches 
$\hat{\Lambda}_1 / \sum_k \hat{\Lambda}_k = 0.339$. This separation 
identifies a coherent structure concentrated at a single frequency. 
An intermediate coherent regime extends up to $f \approx 0.12$, 
where the leading mode retains a meaningful energy advantage over 
the sub-leading modes. Beyond $f \approx 0.14$, the individual 
mode curves converge and the spectrum transitions into an 
incoherent turbulent regime in which energy is distributed nearly 
uniformly across all 19 mode ranks --- the spectral signature of 
broadband, spatially uncorrelated fluctuations.

A localised secondary elevation is clearly visible near $f = 0.12$ 
in both mode 1 and mode 2, breaking the otherwise monotonic decay 
of the leading-mode curve. This bump corresponds to the highest 
rank-separation $\hat{\Lambda}_1 / \hat{\Lambda}_2 = 2.09$ in the 
spectrum (Section~\ref{sec:spectral_summary}) and motivates the 
inclusion of $f = 0.12$ as a target frequency in 
Section~\ref{sec:vertical_cross_sections}.

On the log--log axes, the total energy curve $\sum_k \hat{\Lambda}_k$ 
tracks the $f^{-5/3}$ reference closely from $f \approx 0.05$ up to 
the Nyquist frequency $f = 0.5$. This is broadly consistent with a 
classical Kolmogorov inertial subrange and indicates that the 
small-scale fluctuations remain well-resolved across the full 
sampled bandwidth.

\subsection{Vertical Structure of the Three Leading Modes Across Frequencies}
\label{sec:vertical_cross_sections}

The vertical organisation of the three most energetic SPOD modes 
was examined as a function of frequency. 
Figures~\ref{fig:mode1_xzslices_multifreq}, 
\ref{fig:mode2_xzslices_multifreq}, and~\ref{fig:mode3_xzslices_multifreq} 
present the phase-aligned real parts of $\hat{\Phi}_1$, 
$\hat{\Phi}_2$, and $\hat{\Phi}_3$ at four target frequencies 
($f = 0.02$, $0.06$, $0.12$, and $0.40$), sliced at mid-span 
($y = 4.737$). The four frequencies sample the three spectral 
regimes identified in Section~\ref{sec:spod_spectrum}: the primary 
coherent peak ($f = 0.02$), the intermediate coherent range 
($f = 0.06$ and $f = 0.12$), and the incoherent turbulent reference 
($f = 0.40$).

\subsubsection*{Canonical phase alignment}

SPOD modes are complex-valued, encoding both the spatial amplitude 
and the convective phase of a propagating structure. Because the 
temporal origin of the snapshot ensemble is arbitrary, a direct 
extraction of the real part of $\hat{\Phi}$ would yield an arbitrary 
instantaneous slice of the wave cycle and would obscure direct 
comparison between modes at different frequencies. To visualise 
each mode at a physically reproducible state, a global phase 
rotation $e^{-i\varphi^\star}$ was applied prior to plotting, where 
$\varphi^\star$ is the angle that maximises the $W$-weighted 
real-part energy of the rotated mode. The closed-form solution is
\begin{equation}
    \varphi^\star = \tfrac{1}{2} \arg \left( \sum_n W_n \, \hat{\Phi}_n^2 \right),
    \label{eq:phase_align}
\end{equation}
where $W_n$ are the diagonal entries of the spatial weight matrix 
introduced in Section~\ref{sec:weight_compute} and the sum runs 
over all spatial degrees of freedom. The fraction of total mode 
energy carried by the real part after this alignment is a useful 
diagnostic: a value near $100\%$ corresponds to a standing wave, 
while a value near $50\%$ indicates a purely propagating wave, 
since the energy of a travelling wave is partitioned equally 
between its real and imaginary components. Across all twelve 
mode--frequency combinations considered here, the real-part energy 
fraction lies within $0.3\%$ of the theoretical $50\%$ value, 
confirming that every plotted structure corresponds to a 
downstream-propagating wave.

\begin{figure}[hbt!]
    \centering
    \includegraphics[width=\textwidth]{figures/spod_mode1_multifreq_xzslices.png}
    \caption{Phase-aligned real part of SPOD mode 1 at four 
    frequencies, plotted as vertical (x--z) slices at mid-span 
    ($y = 4.737$). Rows: $U$, $V$, $W$ velocity components. 
    Columns: increasing frequency from left to right. Each panel 
    is normalised to its own maximum amplitude.}
    \label{fig:mode1_xzslices_multifreq}
\end{figure}

\begin{figure}[hbt!]
    \centering
    \includegraphics[width=\textwidth]{figures/spod_mode2_multifreq_xzslices.png}
    \caption{Phase-aligned real part of SPOD mode 2 at four 
    frequencies, plotted as vertical (x--z) slices at mid-span. 
    Layout identical to Figure~\ref{fig:mode1_xzslices_multifreq}.}
    \label{fig:mode2_xzslices_multifreq}
\end{figure}

\begin{figure}[hbt!]
    \centering
    \includegraphics[width=\textwidth]{figures/spod_mode3_multifreq_xzslices.png}
    \caption{Phase-aligned real part of SPOD mode 3 at four 
    frequencies, plotted as vertical (x--z) slices at mid-span. 
    Layout identical to Figure~\ref{fig:mode1_xzslices_multifreq}.}
    \label{fig:mode3_xzslices_multifreq}
\end{figure}

\subsubsection*{Mode 1: dominant coherent structure}

The frequency progression of mode 1 
(Figure~\ref{fig:mode1_xzslices_multifreq}) reveals a systematic 
transformation in spatial character. At $f = 0.02$ the structure 
is large and smooth, with $U$ and $W$ exhibiting paired 
positive/negative regions extending across the full streamwise 
extent of the box. Energy is concentrated in the lower-to-mid 
region, with peak at $z = 0.39$ and centroid $\bar{z} = 0.49$, and 
the component partition $U:V:W = 68.0\%:21.6\%:10.4\%$ confirms 
strong streamwise dominance. At $f = 0.06$ the structure contracts 
to roughly two streamwise periods within $L_x = 3\pi$ and develops 
a forward-leaning inclination visible in the $V$ and $W$ 
components, consistent with mean-shear deformation of coherent 
structures. The mode at $f = 0.12$ exhibits qualitatively different 
behaviour: $U$ weakens substantially while $V$ and especially $W$ 
become concentrated in a thin near-surface layer 
($z \lesssim 0.15$) with high streamwise wavenumber. The component 
partition shifts to $U:V:W = 44.1\%:20.4\%:35.5\%$, with $W$ 
approaching parity with $U$ --- the only frequency in the spectrum 
at which this occurs, and the same frequency that exhibits the 
highest rank-separation $\hat{\Lambda}_1 / \hat{\Lambda}_2 = 2.09$ 
in Section~\ref{sec:spectral_summary}. At $f = 0.40$, all three 
components display fine-scale, spatially disorganised fluctuations 
distributed throughout the wall-normal extent, and the partition 
approaches near-isotropy ($U:V:W \approx 26\%:41\%:33\%$), 
consistent with the incoherent turbulent regime.

\subsubsection*{Modes 2 and 3: independent sub-dominant structures}

At $f = 0.02$, modes 1, 2, and 3 collectively account for 
approximately $65\%$ of the total spectral energy, raising the 
question of whether the sub-dominant modes represent genuinely 
independent structures or are simply variants of mode 1 
(e.g.\ harmonics or phase-shifted copies). The two sub-dominant 
modes (Figures~\ref{fig:mode2_xzslices_multifreq} 
and~\ref{fig:mode3_xzslices_multifreq}) differ from mode 1 in both 
spatial localisation and component partition. Mode 2 places its 
strongest $U$ and $W$ activity in the upper half of the domain 
($z \gtrsim 0.5$), while mode 3 is sharply concentrated in the 
lower quarter ($z \lesssim 0.25$). The component partitions, 
$U:V:W = 72.2\%:17.7\%:10.1\%$ for mode 2 and 
$U:V:W = 76.0\%:14.6\%:9.5\%$ for mode 3, are systematically 
\emph{more} streamwise-dominated than mode 1's 
$68.0\%:21.6\%:10.4\%$. The differing wavelengths and wall-normal 
localisations rule out the harmonic and phase-shifted-copy 
hypotheses, identifying the three modes as distinct coexisting 
structures rather than overtones of a single underlying motion.

At higher frequencies, modes 2 and 3 follow the same general 
trends as mode 1: shorter streamwise wavelengths and broader 
wall-normal distributions in the intermediate range, and 
convergence to a near-identical fine-scale character at 
$f = 0.40$ with statistically indistinguishable component 
partitions ($U:V:W \approx 26\text{--}27\%:40\text{--}41\%:33\%$). 
This convergence reflects the loss of mode-specific spatial 
organisation in the incoherent regime.

\subsubsection*{Component-energy cascade across rank and frequency}

A consistent quantitative trend emerges across the three modes. 
Within the coherent regime, the streamwise share $U$ 
\emph{increases} with mode rank 
($68.0\% \rightarrow 72.2\% \rightarrow 76.0\%$ at $f = 0.02$), 
reflecting an increasingly streak-dominated character in 
higher-rank modes. With increasing frequency, the $U$ share 
\emph{decreases} for every individual mode and the partition 
approaches isotropy. By $f = 0.40$, all three modes share the 
same mildly $V$-biased near-isotropic split. The 
streamwise-anisotropy imprint of the wave-driven shear is 
therefore strongest in the leading mode at the lowest coherent 
frequency and is progressively erased moving toward either higher 
rank or higher frequency.

\subsection{Wall-Normal Energy Distribution and Component-Energy Cascade}
\label{sec:zprofile_uvw}

To complement the qualitative spatial visualisation of the leading 
SPOD modes, two quantitative diagnostics were extracted from the 
mode field. Panel (a) of Figure~\ref{fig:zprofile_uvw} presents 
the wall-normal energy distribution $E(z)$ of mode 1 at the four 
target frequencies, while panel (b) presents the streamwise, 
spanwise, and wall-normal component-energy shares of mode 1 across 
the entire resolved frequency range.

\begin{figure}[hbt!]
    \centering
    \includegraphics[width=\textwidth]{figures/spod_zprofile_uvw.png}
    \caption{(a) Wall-normal energy profile of SPOD mode 1 at the 
    four target frequencies. Profiles are normalised to their own 
    peak and dotted horizontal lines mark the corresponding energy 
    centroids $\bar{z}$. (b) $U$, $V$, and $W$ component-energy 
    shares of mode 1 plotted across all 26 resolved frequency bins. 
    The horizontal grey dotted line marks the equipartition 
    threshold ($1/3$) and the vertical dashed line marks the 
    coherent--incoherent transition identified in 
    Section~\ref{sec:spod_spectrum}.}
    \label{fig:zprofile_uvw}
\end{figure}

\subsubsection*{Wall-normal energy profiles}

The wall-normal distribution of mode 1 varies substantially 
between the four target frequencies, with no two profiles sharing 
the same shape. The relevant integral diagnostics are summarised 
in Table~\ref{tab:zprofile_summary}.

\begin{table}[hbt!]
\centering
\caption{Wall-normal localisation of mode 1 at the four target 
frequencies. $z_{\mathrm{peak}}$ is the height of the energy 
maximum and $\bar{z}$ is the volume-weighted centroid; the 
integral $\int E(z)\,\mathrm{d}z = 1$ in all cases by 
construction, confirming the orthonormality of the modes.}
\label{tab:zprofile_summary}
\begin{tabular}{@{}cccc@{}}
\toprule
\textbf{$f$} & \textbf{$z_{\mathrm{peak}}$} &
\textbf{$\bar{z}$} & \textbf{$\int E(z)\,\mathrm{d}z$} \\
\midrule
0.0200 & 0.390 & 0.488 & 1.000 \\
0.0600 & 0.265 & 0.443 & 1.000 \\
0.1200 & 0.000 & 0.231 & 1.000 \\
0.4000 & 0.247 & 0.386 & 1.000 \\
\bottomrule
\end{tabular}
\end{table}

At $f = 0.02$ the profile is bimodal, with a primary maximum at 
$z_{\mathrm{peak}} = 0.39$ and a near-equal secondary maximum near 
$z \approx 0.95$, yielding a centroid at $\bar{z} = 0.49$. The 
profile at $f = 0.06$ is more conventional, peaking at 
$z_{\mathrm{peak}} = 0.27$ with centroid $\bar{z} = 0.44$, and the 
mode occupies the lower-to-middle portion of the domain. The mode 
at $f = 0.12$ is the only profile peaking at the wave surface 
itself ($z_{\mathrm{peak}} = 0.000$) and drops by an order of 
magnitude within $z \lesssim 0.1$, but a small secondary elevation 
is nevertheless present near the upper boundary; the centroid 
$\bar{z} = 0.23$ confirms the strongly surface-localised nature 
already identified in 
Section~\ref{sec:vertical_cross_sections}. At $f = 0.40$ the 
profile is broader and more diffuse than the coherent profiles, 
with a weak maximum near $z_{\mathrm{peak}} = 0.25$ and a 
secondary elevation in the upper region ($z \approx 0.9$). This 
diffuse distribution is consistent with the spatially incoherent 
nature of the small-scale fluctuations identified in 
Section~\ref{sec:spod_spectrum}. Each frequency therefore samples 
a physically distinct wall-normal organisation rather than a 
monotonic family of progressively detached structures.

\subsubsection*{Component-energy cascade across frequency}

Panel (b) of Figure~\ref{fig:zprofile_uvw} reveals a systematic 
anisotropy cascade in the leading mode. The streamwise share $U$ 
decreases overall from a peak of $69.8\%$ at $f = 0$ to a minimum 
of $26.0\%$ at $f = 0.46$, while the spanwise share $V$ rises from 
$20.5\%$ to a plateau of approximately $40\%$, and the wall-normal 
share $W$ rises from $9.7\%$ to roughly $33\%$. Three distinct 
phases of this cascade can be identified.

In the primary coherent regime ($f \leq 0.04$), the leading mode 
is strongly $U$-dominated, with the partition 
$U:V:W \approx 70\!:\!21\!:\!10$ characteristic of a 
streamwise-elongated streak structure. Across the intermediate 
coherent range ($0.06 \leq f \leq 0.12$), the $U$ share drops 
sharply from $54\%$ to approximately $44\%$ while $V$ and $W$ both 
grow, with $W$ exhibiting an isolated peak of $35.5\%$ at 
$f = 0.12$ --- the only frequency in the spectrum at which $W$ 
exceeds $V$, consistent with the surface-locked, 
wall-normal-dominated structure identified in 
Section~\ref{sec:vertical_cross_sections}. Beyond the 
coherent--incoherent transition at $f \approx 0.14$, all three 
component shares stabilise: the $U$ share crosses the 
equipartition threshold ($1/3$) at $f = 0.20$, indicating the loss 
of streamwise dominance, and from $f \approx 0.30$ onward the 
partition remains approximately 
$U:V:W \approx 27\!:\!40\!:\!33$. The mild $V$-bias in this 
asymptotic state, with $V$ slightly exceeding $W$ and clearly 
exceeding $U$, indicates that the small-scale fluctuations are 
close to but not strictly isotropic, retaining a weak anisotropy 
associated with the wall-bounded geometry.

The cascade summarised here provides a compact quantitative 
description of how energy migrates between velocity components as 
spatial scales decrease: the streamwise-anisotropic imprint of the 
wave-driven shear is strongest at the lowest coherent frequencies 
and is progressively erased as the cascade approaches the 
small-scale turbulent regime.

%==========================================================
\appendix
\section{Complete SPOD Eigenvalue Spectrum}
\label{app:full_spectrum}
%==========================================================

\begin{longtable}{@{}ccccccc@{}}
\caption{Full SPOD spectral data across all resolved frequencies 
($N_f = 26$, $N_{\mathrm{blk}} = 19$).}
\label{tab:full_spod} \\
\toprule
\textbf{idx} & \textbf{$f$} & $\sum \hat{\Lambda}_k$ &
$\hat{\Lambda}_1$ &
$\hat{\Lambda}_1 / \sum \hat{\Lambda}_k$ &
$\hat{\Lambda}_1 / \hat{\Lambda}_2$ &
\textbf{Modes $\geq 90\%$} \\
\midrule
\endfirsthead
\toprule
\textbf{idx} & \textbf{$f$} & $\sum \hat{\Lambda}_k$ &
$\hat{\Lambda}_1$ &
$\hat{\Lambda}_1 / \sum \hat{\Lambda}_k$ &
$\hat{\Lambda}_1 / \hat{\Lambda}_2$ &
\textbf{Modes $\geq 90\%$} \\
\midrule
\endhead
\bottomrule
\endfoot
0  & 0.0000 & 1.6545e+00 & 3.5888e-01 & 0.217 & 1.21 & 8  \\
1  & 0.0200 & 3.4743e+00 & 1.1780e+00 & 0.339 & 1.97 & 10 \\
2  & 0.0400 & 3.0786e+00 & 4.6578e-01 & 0.151 & 1.14 & 14 \\
3  & 0.0600 & 2.2564e+00 & 2.6110e-01 & 0.116 & 1.19 & 16 \\
4  & 0.0800 & 1.7191e+00 & 2.5358e-01 & 0.148 & 1.53 & 16 \\
5  & 0.1000 & 1.3179e+00 & 1.4172e-01 & 0.108 & 1.20 & 16 \\
6  & 0.1200 & 1.1498e+00 & 2.0701e-01 & 0.180 & 2.09 & 16 \\
7  & 0.1400 & 8.7465e-01 & 7.1543e-02 & 0.082 & 1.05 & 16 \\
8  & 0.1600 & 6.9101e-01 & 5.4459e-02 & 0.079 & 1.05 & 17 \\
9  & 0.1800 & 5.7765e-01 & 4.6100e-02 & 0.080 & 1.10 & 17 \\
10 & 0.2000 & 4.9148e-01 & 4.0444e-02 & 0.082 & 1.05 & 17 \\
11 & 0.2200 & 4.1415e-01 & 3.6056e-02 & 0.087 & 1.14 & 17 \\
12 & 0.2400 & 3.5189e-01 & 2.7946e-02 & 0.079 & 1.04 & 17 \\
13 & 0.2600 & 3.0017e-01 & 2.3142e-02 & 0.077 & 1.06 & 17 \\
14 & 0.2800 & 2.5751e-01 & 1.9418e-02 & 0.075 & 1.05 & 17 \\
15 & 0.3000 & 2.2405e-01 & 1.7650e-02 & 0.079 & 1.07 & 17 \\
16 & 0.3200 & 1.9598e-01 & 1.4605e-02 & 0.075 & 1.03 & 17 \\
17 & 0.3400 & 1.7270e-01 & 1.2801e-02 & 0.074 & 1.05 & 17 \\
18 & 0.3600 & 1.5536e-01 & 1.1350e-02 & 0.073 & 1.03 & 17 \\
19 & 0.3800 & 1.3994e-01 & 1.0335e-02 & 0.074 & 1.05 & 17 \\
20 & 0.4000 & 1.2711e-01 & 9.4636e-03 & 0.074 & 1.05 & 17 \\
21 & 0.4200 & 1.1709e-01 & 8.5011e-03 & 0.073 & 1.04 & 17 \\
22 & 0.4400 & 1.1012e-01 & 7.9072e-03 & 0.072 & 1.02 & 17 \\
23 & 0.4600 & 1.0509e-01 & 7.5801e-03 & 0.072 & 1.01 & 17 \\
24 & 0.4800 & 1.0193e-01 & 7.2930e-03 & 0.072 & 1.01 & 17 \\
25 & 0.5000 & 1.0114e-01 & 7.2718e-03 & 0.072 & 1.02 & 17 \\
\midrule
$\Sigma$ & --- & 2.0160e+01 & --- & --- & --- & --- \\
\end{longtable}

\end{document}
