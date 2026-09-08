module EnsembleMCMCCUDAExt

import CUDA
import EnsembleMCMC

EnsembleMCMC._check_kernel_array(initial::CUDA.AnyCuArray) =
    CUDA.functional() || throw(ArgumentError("CUDA is not functional"))

EnsembleMCMC._with_kernel_device(f, initial::CUDA.AnyCuArray) =
    CUDA.device!(f, CUDA.device(initial))

end
