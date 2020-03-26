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
#ifndef COLOR_H
#define COLOR_H

#include <QColor>

#include "util.h"

class Color
{
public:
    Color();
    Color(const QColor &color);
    Color(float L, float a, float b, float alpha = 1.0);

    const QColor getQColor() const { return m_qColor; }

    float getXYZ_X() const { return m_XYZ.X; }
    float getXYZ_Y() const { return m_XYZ.Y; }
    float getXYZ_Z() const { return m_XYZ.Z; }
    void setXYZ(float X, float Y, float Z, float alpha)
    {
        m_XYZ.X = X;
        m_XYZ.Y = Y;
        m_XYZ.Z = Z;
        m_alpha = alpha;
        updateFromXYZ();
    }

    float getxyY_x() const { return m_xyY.x; }
    float getxyY_y() const { return m_xyY.y; }
    float getxyY_Y() const { return m_xyY.Y; }

    float getLab_L() const { return m_Lab.L; }
    float getLab_a() const { return m_Lab.a; }
    float getLab_b() const { return m_Lab.b; }
    void setLab(float L, float a, float b, float alpha)
    {
        m_Lab.L = L;
        m_Lab.a = a;
        m_Lab.b = b;
        m_alpha = alpha;
        updateFromLab();
    }

private:
    // RGB, HSV, HSL, CMYK
    QColor m_qColor;

    float m_alpha;

    // XYZ
    struct
    {
        float X;
        float Y;
        float Z;
    } m_XYZ;

    // xyY
    struct
    {
        float x;
        float y;
        float Y;
    } m_xyY;

    // CIELab
    struct
    {
        float L;
        float a;
        float b;
    } m_Lab;

    void updateFromRGB();
    void updateFromXYZ();
    void updateFromLab();
    void RGBToXYZ();
    void XYZToRGB();
    void XYZToxyY();
    void XYZToLab();
    void LabToXYZ();
};

#endif // COLOR_H
