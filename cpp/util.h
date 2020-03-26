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

#ifndef UTIL_H
#define UTIL_H

#include <QMatrix4x4>

template <typename T>
constexpr bool almostEquals(const T& x, const T& y)
{
    return std::abs(x - y) <= (std::max(std::abs(x), std::abs(y)) *
                               std::numeric_limits<T>::epsilon());
}

template <typename T>
constexpr bool inBounds(const T& v, const T& lo, const T& hi)
{
    return !(v < lo) && (v < hi);
}

constexpr bool inBounds(const QPoint& v, const QPoint& lo, const QPoint& hi)
{
    return v.x() >= lo.x() && v.x() < hi.x() &&
           v.y() >= lo.y() && v.y() < hi.y();
}


// Illuminant = D65 (2°, CIE 1931)
static constexpr QVector3D s_refWhite = QVector3D(0.95047f, 1, 1.08883f);

// RGB-XYZ conversion matrix (sRGB - D65)
static const QMatrix4x4 s_RGBToXYZ = QMatrix4x4(
            0.4124564f,  0.3575761f,  0.1804375f, 0,
            0.2126729f,  0.7151522f,  0.0721750f, 0,
            0.0193339f,  0.1191920f,  0.9503041f, 0,
            0,           0,           0,          1);

// XYZ-RGB conversion matrix (sRGB - D65)
static const QMatrix4x4 s_XYZToRGB = QMatrix4x4(
            3.2404548f, -1.5371389f, -4.9853155f, 0,
           -9.6926639f,  1.8760110f,  4.1556082f, 0,
            5.5643420f, -2.0402585f,  1.0572251f, 0,
            0,           0,           0,          1);

#endif // UTIL_H
