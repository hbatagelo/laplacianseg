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

#ifndef SEEDEDSEGMENTATIONSOFT_H
#define SEEDEDSEGMENTATIONSOFT_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>
#include <QImage>

#include "util.h"

class SeededSegmentationSoft
{
public:
    const QImage* m_image;

    SeededSegmentationSoft(QImage *inputImage,
                           QImage *seedsImage,
                           QColor fgColor = Qt::red,
                           QColor bgColor = Qt::blue,
                           double beta = 500);

    bool Compute();
    void GetOutput(QImage &outputImage, bool asBinary = false);

private:
    // Input image matrix (3 channels)
    Eigen::MatrixXd m_inputImage[3];

    // Sizes of the input image
    int m_pixelCount;
    int m_rows;
    int m_cols;

    double m_beta;

    // Seeds matrix (single value per pixel)
    Eigen::MatrixXd m_seedsMatrix;

    // Sparse matrix storing the weights
    Eigen::SparseMatrix<double> m_sparseWeights;

    // Sparse matrices minimizing the energy functional
    Eigen::SparseMatrix<double> m_sparseD;
    Eigen::SparseMatrix<double> m_sparseIs;
    Eigen::SparseMatrix<double> m_sparseL;

    Eigen::MatrixXd m_solutionMatrix;

    // Convert a matrix into vector
    Eigen::VectorXd Matrix2Vector(Eigen::MatrixXd const &M);

    // Convert a vector into matrix of size rows*cols
    Eigen::MatrixXd Vector2Matrix(Eigen::VectorXd const &V, int rows, int cols);

    // Set up the graph weights and store them in a sparse matrix.
    // Also creates a diagonal sparse matrix with the sum of the
    // weigths in each row.
    void SetUpGraphWeightAndSum();

    // Implement the diagonal sparse matrix with ones, if the pixel belongs
    // to the background or the foreground, and zeros otherwise.
    void CheckIfInSeedsOrNot();

    // Implement the matrix corresponding to sum of the weights
    // minus the weights.
    void DiffBetweenSumAndWeights();

    // Solve the energy functional by minimising using Cholesky
    // factorization algorithm.
    void SolveEnergyFunction();
};

#endif // SEEDEDSEGMENTATION_H
