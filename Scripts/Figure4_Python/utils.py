import numpy as np
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks

def find_thresholds(uPTM, N=10, tol=1e-4, max_iter=100):
    """Find thresholds spanning N sparsity levels."""
    target_sparsities = np.linspace(0, 1, N)
    thresholds = []

    for sparsity in reversed(target_sparsities):
        low, high = np.nanmin(uPTM), np.nanmax(uPTM)
        for _ in range(max_iter):
            mid = (low + high) / 2
            binary_matrix = (uPTM >= mid).astype(int)
            current_sparsity = 1 - np.mean(binary_matrix)

            if abs(current_sparsity - sparsity) < tol:
                break
            elif current_sparsity < sparsity:
                low = mid
            else:
                high = mid
        thresholds.append(mid)

    return np.array(thresholds[::-1])



def run_case(StimulusGains, Inputs, DivNorm1, PositionCorrelation=None, corr_strength=None, num_columns=30):

    # Parameters
    num_neurons = 1000
    Input_Sparse = 0.8
    ConnectStrength = 0.05
    scaleSigma = 0.02
    scaleSigma_position = 3
    num_positions=30

    # Correlation values
    corr_strengths_position_full = np.concatenate(([0], np.linspace(1e-9, 1, 20)))
    corr_strengths_full = np.linspace(0, 1, 21)

    rows_binary = Inputs
    cols_binary = num_columns
    k = int(np.round(num_columns * (1 - Input_Sparse)))

    corr_pos_indices = range(len(corr_strengths_position_full)) if PositionCorrelation is None else [PositionCorrelation]
    corr_indices = range(len(corr_strengths_full)) if corr_strength is None else [corr_strength]

    CPTM_all = [[[None for _ in range(21)] for _ in corr_indices] for _ in corr_pos_indices]
    PTM_all = [[[None for _ in range(21)] for _ in corr_indices] for _ in corr_pos_indices]
    Indices_all = [[[[] for _ in range(21)] for _ in corr_indices] for _ in corr_pos_indices]

    # Position-dependent inputs
    binary_matrix1 = np.zeros((rows_binary, cols_binary), dtype=int)
    for row in range(rows_binary):
        random_indices = np.random.choice(cols_binary, k, replace=False)
        binary_matrix1[row, random_indices] = 1

    x_positions = np.arange(cols_binary)
    MatrixPos = []

    for corr_strength_pos in corr_strengths_position_full:
        sigma = corr_strength_pos * cols_binary / scaleSigma_position
        if corr_strength_pos > 0:
            gaussian_matrix1 = np.zeros((rows_binary, cols_binary))
            for row in range(rows_binary):
                for col in np.where(binary_matrix1[row] == 1)[0]:
                    distances = np.minimum(
                        np.abs(x_positions - col),
                        cols_binary - np.abs(x_positions - col)
                    )
                    gaussian_matrix1[row] += np.exp(-distances**2 / (2 * sigma**2))
            normalized = gaussian_matrix1 / np.nanmax(gaussian_matrix1)
            normalized *= StimulusGains
        else:
            normalized = binary_matrix1 * StimulusGains
        MatrixPos.append(normalized)

    # Random input connectivity
    base_binary_matrix = (np.random.rand(num_neurons, Inputs) < ConnectStrength).astype(float)
    uniform_random_values = np.random.uniform(0, 1, base_binary_matrix.shape)

    for Q_idx, Q in enumerate(corr_pos_indices):
        Matrix = MatrixPos[Q]

        for j_idx, j in enumerate(corr_indices):
            corr_strength_val = corr_strengths_full[j]

            if corr_strength_val == 0:
                correlated_binary_matrix = base_binary_matrix.copy()
            else:
                smoothed_matrix = gaussian_filter(
                    base_binary_matrix, sigma=(corr_strength_val / scaleSigma, 0), mode="constant"
                )
                threshold = np.percentile(smoothed_matrix, 100 * (1 - ConnectStrength))
                correlated_binary_matrix = (smoothed_matrix >= threshold).astype(float)

            correlated_binary_matrix = np.where(
                correlated_binary_matrix == 1, uniform_random_values, 0
            )

            preuPTM = np.dot(correlated_binary_matrix, Matrix)

            # Divisive normalization
            if DivNorm1 == 1:
                uPTM_Max = preuPTM / np.nanmax(preuPTM, axis=0)
            else:
                uPTM_Max = preuPTM / np.nanmax(preuPTM)

            thresholds_max = find_thresholds(uPTM_Max, N=21)

            for argmax_idx, thresh in enumerate(thresholds_max):
                df = pd.DataFrame(uPTM_Max)
                preferred_columnsNew = np.array(df.idxmax(axis=1))
                Indices_all[Q_idx][j_idx][argmax_idx] = preferred_columnsNew

                tPTM = np.where(uPTM_Max >= thresh, uPTM_Max, 0)
                SparsityM = (uPTM_Max >= thresh).astype(int)
                CPTM_all[Q_idx][j_idx][argmax_idx] = tPTM
                PTM_all[Q_idx][j_idx][argmax_idx] = SparsityM

    return CPTM_all, PTM_all, Indices_all



def simulate_neural_activity(MultiMatrix_Arena, MultiMatrix_Room, WeightArena, WeightRoom, T, velocity, Sigma, phase_shift_rate, DivNorm2, Thresh, obj_placements):
    """Simulate combined arena- and room-frame activity."""
    N, P = MultiMatrix_Arena.shape

    positions_arena = np.floor(np.cumsum(np.full(T, velocity)) % P).astype(int)

    # Room frame shifts relative to the arena frame
    positions_room = ((positions_arena - np.floor(np.arange(T) * phase_shift_rate) % P) % P).astype(int)

    activations = np.zeros((N, T))

    for t in range(T):
        pos_arena = positions_arena[t]
        pos_room = positions_room[t]

        activity = WeightArena*MultiMatrix_Arena[:, pos_arena].copy()
        activity += WeightRoom * MultiMatrix_Room[:, pos_room]

        # Normalize and sparsify combined activity
        if DivNorm2 == 1:
            activity /= np.max(activity, axis=0)
        else:
            activity /= np.max(activity)

        thresh = find_thresholds(activity, N=21)
        activity = np.where(activity >= thresh[Thresh], activity, 0)

        activity += np.random.normal(0, Sigma, size=N)
        activity = np.clip(activity, 0, None)

        activations[:, t] = activity

    return activations, positions_arena, positions_room


def compute_tuning_curve(activations, positions, P):
    N, T = activations.shape
    tuning_curves = np.zeros((N, P))

    for p in range(P):
        mask = positions == p
        if np.any(mask):
            tuning_curves[:, p] = np.mean(activations[:, mask], axis=1)

    return tuning_curves

def sort_neurons_by_peak(tuning_curves):
    peak_positions = np.argmax(tuning_curves, axis=1)
    sorted_indices = np.argsort(peak_positions)

    return sorted_indices, peak_positions

def compute_rayleigh_vector_length(tuning_curves, P):
    N, _ = tuning_curves.shape
    angles = np.linspace(0, 2 * np.pi, P, endpoint=False)

    x_component = np.sum(tuning_curves * np.cos(angles), axis=1)
    y_component = np.sum(tuning_curves * np.sin(angles), axis=1)

    r = np.sqrt(x_component**2 + y_component**2) / np.sum(tuning_curves, axis=1)

    return r

def compute_skaggs_information(tuning_curves, P):
    N, _ = tuning_curves.shape
    p_x = 1 / P

    mean_firing_rates = np.mean(tuning_curves, axis=1, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        info_content = np.nansum(
            p_x * (tuning_curves / mean_firing_rates) * np.log2(tuning_curves / mean_firing_rates),
            axis=1
        )

    info_content[np.isnan(info_content)] = 0
    return info_content

def circular_difference(pos_A, pos_C, P):
    diff = (pos_C - pos_A + P // 2) % P - P // 2
    return diff

def circular_center_of_mass(tuning_curves):
    """Circular center of mass for each tuning curve."""
    num_positions = tuning_curves.shape[1]
    angles = np.linspace(0, 2 * np.pi, num_positions, endpoint=False)

    x = np.sum(tuning_curves * np.cos(angles), axis=1)
    y = np.sum(tuning_curves * np.sin(angles), axis=1)

    com = np.arctan2(y, x)
    com = (com + 2 * np.pi) % (2 * np.pi)

    com = com * (num_positions / (2 * np.pi))

    return com