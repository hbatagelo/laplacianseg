/****************************************************************************
**
** Copyright (C) 2020 Harlen Batagelo <hbatagelo@gmail.com> and
**                    João Paulo Gois <jpgois@gmail.com>.
**
** This file is part of the implementation of the paper
** 'Laplacian Coordinates: Theory and Methods for Seeded Image Segmentation'
** by Wallace Casaca, João Paulo Gois, Harlen Batagelo, Gabriel Taubin and
** Luis Gustavo Nonato.
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program. If not, see <http://www.gnu.org/licenses/>.
**
****************************************************************************/

#ifndef SEEDEDSPSEGMENTATION_H
#define SEEDEDSPSEGMENTATION_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>
#include <QImage>

#include "slic.h"
#include "opencv2/opencv.hpp"

class SeededSPSegmentation
{
public:
    SeededSPSegmentation(QImage *inputImage,
                         QImage *seedsImage,
                         int superpixelSize,
                         double compactness);

    bool Compute(bool useHardSeeds = false,
                 QColor fgColor = Qt::red,
                 QColor bgColor = Qt::blue,
                 double betaF = 88,
                 double betaD = -1);

    void GetOutput(QImage &outputImage, bool asBinary = false);

private:
    QImage *m_image; // Input image
    QImage *m_seeds; // Input seeds image

    SLIC m_slic;

    // Seeds vector (single value per superpixel)
    Eigen::VectorXd m_seedsVector;

    // Sparse matrix storing the weights
    Eigen::SparseMatrix<double> m_sparseWeights;

    // All sparse matrices minimizing the energy functional
    Eigen::SparseMatrix<double> m_sparseD;
    Eigen::SparseMatrix<double> m_sparseIs;
    Eigen::SparseMatrix<double> m_sparseL;

    Eigen::VectorXd m_solutionVector;

    bool m_useHardSeeds;

    // Returns the beta which results in the largest variation of weights
    double ComputeDistanceBeta(double betaMin = 0.1,
                               double betaMax = 0.5,
                               double betaStep = 0.01,
                               bool global = false);

    // Set up the graph weights and store them in a sparse matrix
    // and deliver also a diagonal Sparse Matrix with the sum of the
    // weigths in each row
    void SetUpGraphWeightAndSum(double betaFeatures, double betaDistance);

    // Implement the diagonal sparse matrix with ones, if the pixel belongs
    // to the background or the foreground, and zeros otherwise.
    void CheckIfInSeedsOrNot();

    // Implement the matrix corresponding to sum of the weights
    // minus the weights
    void DiffBetweenSumAndWeights();

    // Solve the energy functional by minimising using Cholesky
    // factorization algorithm
    void SolveSoftLinearSystem();
    void SolveHardLinearSystem();

    bool GetFirstNonZeroPixel(cv::Mat &image, cv::Mat &segImage, cv::Point &point);
};

#endif // SEEDEDSPSEGMENTATION_H
