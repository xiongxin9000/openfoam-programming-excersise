/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    filedmanipulation

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
 
// 声明一个函数，后面有该函数的具体实现过程
// 输入时间t，空间坐标x，参考点x0以及缩放因子scale
scalar calculatedPressure(scalar t, vector x, vector x0, scalar scale);
 
int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
 
    Info << "read transportProperties" << endl;
    // 从transportProperties文件中读取数据
    // 先定义一个IOdictionary对象，其构造函数参数为一个IOobject对象
    IOdictionary transportProperties(
        IOobject(
            "transportProperties",           // 字典文件名
            runTime.constant(),              // 字典文件所在路径，这里为constant文件夹下
            mesh,                            // 一个objectRegistry类对象，这里没什么用
            IOobject::MUST_READ_IF_MODIFIED, // 如果文件被修改，则必须重新读取
            IOobject::NO_WRITE               //表示文件为只读
            ));
 
    // 定义一个dimensionedScalar变量nu
    // nu有量纲，其量纲dimViscosity等同于(0,2,-1,0,0,0,0),单位为m2/s
    dimensionedScalar nu(
        "nu",         //指定名称
        dimViscosity, //指定scalar的量纲
        // 旧版本或org版本的写法：transportProperties.lookup("nu")
        transportProperties // com新版本写法，自动根据名称在字典中搜索
    );
    Info << "nu is :" << nu << endl;
 
    //-------------从p文件中读取压力场p--------
    // 不需要查找关键字，因为整个文件都是关于压力场的数据
    Info << "read pressure" << endl;
    // 定义一个标量场p，无需指定量纲，因为其量纲已经在相应的文件中指定了
    volScalarField p(
        IOobject(
            "p",                //指定名称
            runTime.timeName(), // 获取当前时间
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh);
 
    // -------从U文件中读取速度场---------
    Info << "read velocity" << endl;
    // 定义一个向量场U
    volVectorField U(
        IOobject(
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh);
 
    // 定义一个场向量
    const vector originVector(0.05, 0.05, 0.005);
 
    // 计算向量originVector与各网格中心的距离中的最大值
    // 利用dimensionedVector将vector转换为具有长度量纲的向量
    // mag函数计算向量的模
    // value函数将数值转换为无量纲标量值
    const scalar rFarCell = max(mag(dimensionedVector("x0",dimLength,originVector)-mesh.C())).value();
     
    // 在案例迭代过程中进行场变量计算
    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        // 在所有的网格上循环
        for (label cellI = 0; cellI < mesh.C().size(); cellI++)
        {
            // 根据自定义的函数计算压力值
            p[cellI] = calculatedPressure(runTime.time().value(), mesh.C()[cellI],
                                          originVector, rFarCell);
        }
 
        // 计算压力梯度，并转换为速度，由于量纲不一致，所以需要乘以一个时间量纲进行转换
        // 注意在不可压缩求解器中，压力的量纲为[0 2 -2 0 0 0 0]，压力梯度量纲[0 1 -2 0 0 0 0]
        // 因此压力梯度乘上时间得到速度量纲[0 1 -1 0 0 0 0]
        U = fvc::grad(p) * dimensionedScalar("tmp", dimTime, 1.0);
 
        // 将数据写入到文件中，标记为AUTO_WRITE的场数据才可被写入
        runTime.write();
    }
 
    Info << "simulation end" << endl;
 
    return 0;
}
 
scalar calculatedPressure(scalar t, vector x, vector x0, scalar scale)
{
    scalar r(mag(x - x0) / scale);
    // 分母加1e-12是为了防止除以零
    scalar rR(1/(r+1e-12));
    scalar f(1.0);
    // 返回一个以x0为中心的正弦分布的压力
    // 圆的半径随着时间逐渐增大
    return f * Foam::cos(2 * Foam::constant::mathematical::pi * rR * t);
}
