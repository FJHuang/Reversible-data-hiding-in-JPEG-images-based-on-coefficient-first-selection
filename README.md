# Reversible-data-hiding-in-JPEG-images-based-on-coefficient-first-selection

DOI: 10.1016/j.sigpro.2022.108639

## Abstract

Joint photographic experts group (JPEG) is the most commonly used image format on the Internet, and using reversible data hiding (RDH) technology to verify its authenticity and integrity is of great significance. So far, many RDH algorithms for JPEG image have been proposed, which improve the performance of the algorithm via block sorting, frequency selection and other adaptive strategies. In this paper, a new JPEG RDH method based on coefficient-first selection strategy is proposed. In our new method, the distribution models of alternating current (AC) coefficients belonging to different discrete cosine transform (DCT) coefficient blocks are constructed first, and then the distortion value of each AC coefficient can be calculated through the estimated probability density function of the constructed coefficient model. In the embedding process, the coefficients with low distortion values will be preferentially selected for data hiding. Various experimental results demonstrate that compared with the state-of-the-art JPEG RDH methods, our proposed may have better visual quality and less increase of file storage size.

## Code

optimal_1D.m: 1D histogram shifting based algorithms

optimal_2D.m: 2D histogram shifting based algorithms
