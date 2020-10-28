{
    word alphaScheme("div(phi,alpha)"); //fcSchemes 中的散度格式
    word alpharScheme("div(phirb,alpha)");//Gauss vanLeer   和  Gauss linear  MULES就是FCT的变种，区别在于限制器的求解和原始文献有点不同。使用FCT的时候，需要指定一种高阶格式，比如linear或者vanLeer

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    {
//智能指针        
        tmp<fv::ddtScheme<scalar>> tddtAlpha  //tmp 类是 OpenFOAM 中用来封装对象的一个类
        (
            fv::ddtScheme<scalar>::New   //
            (
                mesh,
                mesh.ddtScheme("ddt(alpha)")
            )
        );
        const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();  
	//ddtAlpha是一个初始化为tddtAlpha()的fv::ddtScheme<scalar>型引用

        if			
//判断ddtAlpha的格式 EulerDdtScheme localEulerDdtScheme CrankNicolsonDdtScheme
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
         || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
        )
        {
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
        {
            if (nAlphaSubCycles > 1)
            {
                FatalErrorInFunction
                    << "Sub-cycling is not supported "
                       "with the CrankNicolson ddt scheme"
                    << exit(FatalError);
            }

            if
            (
                alphaRestart
             || mesh.time().timeIndex() > mesh.time().startTimeIndex() + 1
            )
            {
                ocCoeff =
                    refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
                   .ocCoeff();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    // Set the time blending factor, 1 for Euler 设置均化系数，euler格式下为1
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);
    // 在fvSolution中   nAlphaCorr      1;   nAlphaSubCycles 3;      cAlpha          1;
    // Standard face-flux compression coefficient
    surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));  
    //magSf 网格面的大小  mag 取模
    //压缩项的速度 压缩因子cAlpha=1定义在fvSolution中

    // Add the optional isotropic compression contribution 增加其他的各项同性压缩贡献因子
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    // Add the optional shear compression contribution 增加其他的剪切压缩贡献因子
    if (scAlpha > 0)
    {
        phic +=
            scAlpha*mag(mesh.delta() & fvc::interpolate(symm(fvc::grad(U))));
    }


    surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.) 所有边界出的通量为0
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }
//智能指针phiCN
    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux计算偏离体心通量
    // scalar cnCoeff = 1.0/(1.0 + ocCoeff);
    if (ocCoeff > 0)
    {
        phiCN = cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime();
    }
//fvSolution 中 MULESCorr   yes;求出 下一时间步的alpha1
    if (MULESCorr)   
    {
        //定义zeroField Su;zeroField Sp;zeroField divU;"alphaSuSp.H" 定义一个为0 的场
        #include "alphaSuSp.H"   //申明源项，全部为zeroField类型

        fvScalarMatrix alpha1Eqn
        (
            (              //? 表示条件运算符 if (LTS) {fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)  else  fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)}
//  其中，Exp1、Exp2 和 Exp3 是表达式。请注意冒号的使用和位置。? : 表达式的值取决于 Exp1 的计算结果。如果 Exp1 为真，则计算 Exp2 的值，且 Exp2 的计算结果则为整个 ? : 表达式的值。如果 Exp1 为假，则计算 Exp3 的值，且 Exp3 的计算结果则为整个 ? : 表达式的值。          
                LTS             
              ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
              : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiCN,
                upwind<scalar>(mesh, phiCN)
            ).fvmDiv(phiCN, alpha1)
       // - fvm::Sp(fvc::ddt(dimensionedScalar(dimless, 1), mesh)
       //           + fvc::div(phiCN), alpha1)
         ==                             //检查两个操作数的值是否相等，如果相等，则条件为真
            Su + fvm::Sp(Sp + divU, alpha1)
        );

        alpha1Eqn.solve();    //在fvMatric.C  fvMatrixSolve.C中 用不同的求解线性矩阵的算法 和不同的惨差求解

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
//智能指针phiCN
        tmp<surfaceScalarField> talphaPhi1UD(alpha1Eqn.flux());
        alphaPhi10 = talphaPhi1UD();

        if (alphaApplyPrevCorr && talphaPhi1Corr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct        //correct调用 
            (
                geometricOneField(),   
                alpha1,
                alphaPhi10,
                talphaPhi1Corr0.ref(),
                oneField(),
                zeroField()
            );

            alphaPhi10 += talphaPhi1Corr0();
        }

        // Cache the upwind-flux
        talphaPhi1Corr0 = talphaPhi1UD;

        alpha2 = 1.0 - alpha1;

        mixture.correct();   //相变修正
    }
/*

template
<
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void Foam::MULES::correct
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    const fvMesh& mesh = psi.mesh();

    if (fv::localEulerDdt::enabled(mesh))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);

        limitCorr  //先调用limiterCorr计算 lamda，记录在类型为scalarField的变量allLambda中，
        (
            rDeltaT,
            rho,
            psi,
            phi,
            phiCorr,
            Sp,
            Su,
            psiMax,
            psiMin
        );

        correct(rDeltaT, rho, psi, phiCorr, Sp, Su);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();

        limitCorr
        (
            rDeltaT,
            rho,
            psi,
            phi,
            phiCorr,
            Sp,
            Su,
            psiMax,
            psiMin
        );

        correct(rDeltaT, rho, psi, phiCorr, Sp, Su);
    }
}


*/

// fvSolution 中  nAlphaCorr = 2  ;
    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)   
    {
        //定义zeroField Su;zeroField Sp;zeroField divU;"alphaSuSp.H" 定义一个为0 的场
        #include "alphaSuSp.H"    // 声明源项，全部为zeroField类型

        surfaceScalarField phir(phic*mixture.nHatf());   //nHatf为自由表面的单位法向量
//智能指针
        tmp<surfaceScalarField> talphaPhi1Un              //组建alpha通量，采用时间项若采用Euler，则cnCoeff= 1.0
        (
            fvc::flux   //压缩项
            (
                phiCN(),
                cnCoeff*alpha1 + (1.0 - cnCoeff)*alpha1.oldTime(),
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        if (MULESCorr)                  //MULESCorr
        {
            tmp<surfaceScalarField> talphaPhi1Corr(talphaPhi1Un() - alphaPhi10);
            volScalarField alpha10("alpha10", alpha1);

            MULES::correct              //MULES
            (
                geometricOneField(),
                alpha1,
                talphaPhi1Un(),
                talphaPhi1Corr.ref(),
                Sp,
                (-Sp*alpha1)(),
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi10 += talphaPhi1Corr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                alphaPhi10 += 0.5*talphaPhi1Corr();
            }
        }
        else
        {
            alphaPhi10 = talphaPhi1Un;

            MULES::explicitSolve   //相输运方程的压缩项 MULES
            (
                geometricOneField(),
                alpha1,
                phiCN,
                alphaPhi10,
                Sp,
                (Su + divU*min(alpha1(), scalar(1)))(),
                oneField(),
                zeroField()
            );
        }

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        talphaPhi1Corr0 = alphaPhi10 - talphaPhi1Corr0;
        talphaPhi1Corr0.ref().rename("alphaPhi1Corr0");
    }
    else
    {
        talphaPhi1Corr0.clear();
    }

    #include "rhofs.H"

    if
    (
        word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::EulerDdtScheme<vector>::typeName
     || word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::localEulerDdtScheme<vector>::typeName
    )
    {
        rhoPhi = alphaPhi10*(rho1f - rho2f) + phiCN*rho2f;
    }
    else
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhi10 =
                (alphaPhi10 - (1.0 - cnCoeff)*alphaPhi10.oldTime())/cnCoeff;
        }

        // Calculate the end-of-time-step mass flux
        rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f;
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
