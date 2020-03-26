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

#include <cmath>

#include "color.h"
#include "util.h"

Color::Color()
{
    m_alpha = 1.0;
    m_qColor = QColor(255, 255, 255, 255);
    updateFromRGB();
}

Color::Color(const QColor &color) : m_qColor(color)
{
    m_alpha = static_cast<float>(color.alphaF());
    updateFromRGB();
}

Color::Color(float L, float a, float b, float alpha)
{
    setLab(L, a, b, alpha);
}

void Color::updateFromRGB()
{
    RGBToXYZ();
    XYZToxyY();
    XYZToLab();
}

void Color::updateFromXYZ()
{
    XYZToRGB();
    XYZToxyY();
    XYZToLab();
}

void Color::updateFromLab()
{
    LabToXYZ();
    XYZToxyY();
    XYZToRGB();
}

void Color::RGBToXYZ()
{
    QVector3D v(static_cast<float>(m_qColor.redF()),
                static_cast<float>(m_qColor.greenF()),
                static_cast<float>(m_qColor.blueF()));

    // Inverse sRGB companding
    if (v.x() > 0.04045f)
        v.setX(std::pow((v.x() + 0.055f) / 1.055f, 2.4f));
    else
        v.setX(v.x() / 12.92f);
    if (v.y() > 0.04045f)
        v.setY(std::pow((v.y() + 0.055f) / 1.055f, 2.4f));
    else
        v.setY(v.y() / 12.92f);
    if (v.z() > 0.04045f)
        v.setZ(std::pow((v.z() + 0.055f) / 1.055f, 2.4f));
    else
        v.setZ(v.z() / 12.92f);

    v = s_RGBToXYZ.mapVector(v);
    m_XYZ.X = v.x();
    m_XYZ.Y = v.y();
    m_XYZ.Z = v.z();
}

void Color::XYZToRGB()
{
    QVector3D v(m_XYZ.X, m_XYZ.Y, m_XYZ.Z);

    v = s_XYZToRGB.mapVector(v);

    // sRGB companding
    if (v.x() > 0.0031308f)
        v.setX(1.055f * (std::pow(v.x(), 1 / 2.4f)) - 0.055f);
    else
        v.setX(12.92f * v.x());
    if (v.y() > 0.0031308f)
        v.setY(1.055f * (std::pow(v.y(), 1 / 2.4f)) - 0.055f);
    else
        v.setY(12.92f * v.y());
    if (v.z() > 0.0031308f)
        v.setZ(1.055f * (std::pow(v.z(), 1 / 2.4f)) - 0.055f);
    else
        v.setZ(12.92f * v.z());

    m_qColor.setRedF(static_cast<qreal>(std::clamp(v.x(), 0.0f, 1.0f)));
    m_qColor.setGreenF(static_cast<qreal>(std::clamp(v.y(), 0.0f, 1.0f)));
    m_qColor.setBlueF(static_cast<qreal>(std::clamp(v.z(), 0.0f, 1.0f)));
    m_qColor.setAlphaF(static_cast<qreal>(m_alpha));
}

void Color::XYZToxyY()
{
    float den = m_XYZ.X + m_XYZ.Y + m_XYZ.Z;
    if (!almostEquals(den, 0.0f))
    {
        m_xyY.x = m_XYZ.X / den;
        m_xyY.y = m_XYZ.Y / den;
    }
    else
    {
        m_xyY.x = s_refWhite.x();
        m_xyY.y = s_refWhite.y();
    }
    m_xyY.Y = m_XYZ.Y;
}

void Color::XYZToLab()
{
    QVector3D f(m_XYZ.X / s_refWhite.x(),
                m_XYZ.Y / s_refWhite.y(),
                m_XYZ.Z / s_refWhite.z());

    if (f.x() > 0.008856f)
        f.setX(std::pow(f.x(), 1.0f / 3.0f));
    else
        f.setX((7.787f * f.x()) + (16.0f / 116.0f));

    if (f.y() > 0.008856f)
        f.setY(std::pow(f.y(), 1.0f / 3.0f));
    else
        f.setY((7.787f * f.y()) + (16.0f / 116.0f));

    if (f.z() > 0.008856f)
        f.setZ(std::pow(f.z(), 1.0f / 3.0f));
    else
        f.setZ((7.787f * f.z()) + (16.0f / 116.0f));

    m_Lab.L = (116.0f * f.y() - 16.0f);
    m_Lab.a = 500.0f * (f.x() - f.y());
    m_Lab.b = 200.0f * (f.y() - f.z());
}

void Color::LabToXYZ()
{
    m_XYZ.Y = (m_Lab.L + 16.0f) / 116.0f;
    m_XYZ.X = m_Lab.a / 500.0f + m_XYZ.Y;
    m_XYZ.Z = m_XYZ.Y - m_Lab.b / 200.0f;

    if (std::pow(m_XYZ.Y, 3) > 0.008856)
        m_XYZ.Y = std::pow(m_XYZ.Y, 3.0f);
    else
        m_XYZ.Y = (m_XYZ.Y - 16.0f / 116.0f) / 7.787f;

    if (std::pow(m_XYZ.X, 3) > 0.008856)
        m_XYZ.X = std::pow(m_XYZ.X, 3.0f);
    else
        m_XYZ.X = (m_XYZ.X - 16.0f / 116.0f) / 7.787f;

    if (std::pow(m_XYZ.Z, 3) > 0.008856)
        m_XYZ.Z = std::pow(m_XYZ.Z, 3.0f);
    else
        m_XYZ.Z = (m_XYZ.Z - 16.0f / 116.0f) / 7.787f;

    m_XYZ.X = s_refWhite.x() * m_XYZ.X;
    m_XYZ.Y = s_refWhite.y() * m_XYZ.Y;
    m_XYZ.Z = s_refWhite.z() * m_XYZ.Z;
}
