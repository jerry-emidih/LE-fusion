# LE-fusion
MATLAB code for graph Laplacian-based data fusion tasks

We use Laplacian eigenmaps (LE) for data representation, based on the original paper here https://ieeexplore.ieee.org/abstract/document/6789755.
I will be updating this with an example script soon...

The general idea is that data points which are expressed in two modalities can be used to learn mappings between the modalities. This is the continuance of work from Alexander Cloninger, Wojciech Czaja, and Timothy Doster <https://ieeexplore.ieee.org/document/6946659>, with their code (with minimal, documented changes) given in the oldCode folder. In particular, this relies on their original pre-image algorithm for LE, described in this paper <https://iopscience.iop.org/article/10.1088/1361-6420/aa5489/pdf>.

# BW to Color Example
For example, if we want to transfer black and white (BW) pixels in an image to color pixels, then we can create LE embeddings for the two sets of pixel vectors. This assumes that we have some pixels for which we know the BW and color values. Then, we create a mapping from the BW embedding to the color embedding. This can be done using a rotation operator described by Coifman and Hirn <https://www.sciencedirect.com/science/article/pii/S1063520313000225>, with code in the oldCode folder. This mapping is then applied to new, out-of-sample BW points, sending them into the color LE space. From here, they can be lifted into RGB color space using the pre-image of the LE color embedding. This provides a modality transfer of BW pixels of an image to their approximate color values within the context of that same image.  
