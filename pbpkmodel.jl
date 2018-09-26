using PuMaS

pbpkmodel = @model begin
    @param begin
        GER ∈ ConstDomain(0.066)
        ρ ∈ ConstDomain(5e-6)
        r ∈ ConstDomain(1)
        T ∈ ConstDomain(3e-5)
        D ∈ ConstDomain(1e-4)
        SST ∈ ConstDomain(5.5)
        kilST ∈ ConstDomain(0.5)
        kaST ∈ ConstDomain(14040.00076)
        kaGU ∈ ConstDomain(14040.000063)
        kt ∈ ConstDomain(0.035)
        SGU1 ∈ ConstDomain(5.5)
        SGU2 ∈ ConstDomain(5.5)
        SGU3 ∈ ConstDomain(5.5)
        SGU4 ∈ ConstDomain(5.5)
        SGU5 ∈ ConstDomain(5.5)
        SGU6 ∈ ConstDomain(5.5)
        SGU7 ∈ ConstDomain(5.5)
        kilGU1 ∈ ConstDomain(0.5 )
        kilGU2 ∈ ConstDomain(0.5)
        kilGU3 ∈ ConstDomain(0.5)
        kilGU4 ∈ ConstDomain(0.5)
        kilGU5 ∈ ConstDomain(0.5)
        kilGU6 ∈ ConstDomain(0.5)
        kilGU7 ∈ ConstDomain(0.5)
        EHR ∈ ConstDomain(0 )
        kbil ∈ ConstDomain(0.0)
        VLI ∈ ConstDomain(1690)
        Kp ∈ ConstDomain(1.3)
        ktCO ∈ ConstDomain(0.0007)
        SCO ∈ ConstDomain(5.5)
        VCO ∈ ConstDomain(700)
        kilCO ∈ ConstDomain(0.0007)
        kaCO ∈ ConstDomain(14040.0000542)
        CP ∈ ConstDomain(0)
        QLU ∈ ConstDomain(5233)
        VLU ∈ ConstDomain(1172)
        VST1 ∈ ConstDomain(50)
        VST2 ∈ ConstDomain(154)
        VGU ∈ ConstDomain(1650)
        VAR ∈ ConstDomain(1698)
        AIR ∈ ConstDomain(0.0)
        VVE ∈ ConstDomain(3396)
        VIR ∈ ConstDomain(0.0)
        QBR ∈ ConstDomain(700)
        VBR ∈ ConstDomain(1450)
        QLI ∈ ConstDomain(1650)
        QKI ∈ ConstDomain(1100)
        QHR ∈ ConstDomain(150)
        VHR ∈ ConstDomain(310)
        QMU ∈ ConstDomain(750)
        VMU ∈ ConstDomain(35000)
        QAD ∈ ConstDomain(260)
        VAD ∈ ConstDomain(10000)
        QSK ∈ ConstDomain(300)
        VSK ∈ ConstDomain(7800)
        QBO ∈ ConstDomain(250)
        VBO ∈ ConstDomain(4579)
        QTH ∈ ConstDomain(80)
        VTH ∈ ConstDomain(29)
        QST ∈ ConstDomain(38)
        QGU ∈ ConstDomain(1100)
        Ker ∈ ConstDomain(0.0)
        QPA ∈ ConstDomain(133)
        VPA ∈ ConstDomain(77)
        QSP ∈ ConstDomain(77)
        VSP ∈ ConstDomain(192)
        CLint ∈ ConstDomain(0)
        QHA ∈ ConstDomain(302)
        VKI ∈ ConstDomain(280)
        R ∈ ConstDomain(1)
    end

    @random begin
        η ~ MvNormal(eye(2))
    end

    @dynamics begin
        #Absorption compartments
        #Stomach compartment
        AUNDST' = -GER * AUNDST - (((3*D)/(ρ*r*T)) * AUNDST* (SST - (ADIST/VST1)))
        ADIST' = -GER * ADIST + (((3*D)/(ρ*r*T)) *AUNDST*(SST - (ADIST/VST1))) - kilST* ADIST -kaST *ADIST
        ADEGST' = -GER * ADEGST + kilST * ADIST
        AABSST' = kaST * ADIST

        #GU1 small intestinal compartment
        AUNDGU1' = GER * AUNDST - kt * AUNDGU1 - ((3*D)/(ρ*r*T)) * AUNDGU1 *(SGU1 - (ADISGU1/VGU))
        ADISGU1' = GER * ADIST - kt * ADISGU1 + ((3*D)/(ρ*r*T)) * AUNDGU1 *(SGU1 - (ADISGU1/VGU)) - kilGU1*ADISGU1 - kaGU*ADISGU1 + (EHR*kbil*CLI *VLI)/( Kp)
        ADEGGU1' = GER * ADEGST - kt * ADEGGU1 + kilGU1 * ADISGU1
        AABSGU1' = kaGU * ADISGU1

        # Other small intestinal compartments (GU2–GU7)
        AUNDGU2' = kt * AUNDGU1 - kt * AUNDGU2 -((3*D)/(ρ*r*T)) * AUNDGU2 *(SGU2 - (ADISGU2/VGU))
        ADISGU2' = kt * ADISGU1 - kt * ADISGU2 +((3*D)/(ρ*r*T)) * AUNDGU2 *(SGU2 - (ADISGU2/VGU)) - kilGU2*ADISGU2 - kaGU*ADISGU2
        ADEGGU2' = kt*ADEGGU1 - kt*ADEGGU2 + kilGU2 * ADISGU2
        AABSGU2' = kaGU * ADISGU2

        AUNDGU3' = kt * AUNDGU2 - kt * AUNDGU3 -((3*D)/(ρ*r*T)) * AUNDGU3 *(SGU3 - (ADISGU3/VGU))
        ADISGU3' = kt * ADISGU2 - kt * ADISGU3 +((3*D)/(ρ*r*T)) * AUNDGU3 *(SGU3 - (ADISGU3/VGU)) - kilGU3*ADISGU3 - kaGU*ADISGU3
        ADEGGU3' = kt*ADEGGU2 - kt*ADEGGU3 + kilGU3 * ADISGU3
        AABSGU3' = kaGU * ADISGU3

        AUNDGU4' = kt * AUNDGU3 - kt * AUNDGU4 -((3*D)/(ρ*r*T)) * AUNDGU4 *(SGU4 - (ADISGU4/VGU))
        ADISGU4' = kt * ADISGU3 - kt * ADISGU4 +((3*D)/(ρ*r*T)) * AUNDGU4 *(SGU4 - (ADISGU4/VGU)) - kilGU4*ADISGU4 - kaGU*ADISGU4
        ADEGGU4' = kt*ADEGGU3 - kt*ADEGGU4 + kilGU4 * ADISGU4
        AABSGU4' = kaGU * ADISGU4

        AUNDGU5' = kt * AUNDGU4 - kt * AUNDGU5 -((3*D)/(ρ*r*T)) * AUNDGU5 *(SGU5 - (ADISGU5/VGU))
        ADISGU5' = kt * ADISGU4 - kt * ADISGU5 +((3*D)/(ρ*r*T)) * AUNDGU5 *(SGU5 - (ADISGU5/VGU)) - kilGU5*ADISGU5 - kaGU*ADISGU5
        ADEGGU5' = kt*ADEGGU4 - kt*ADEGGU5 + kilGU5 * ADISGU5
        AABSGU5' = kaGU * ADISGU5

        AUNDGU6' = kt * AUNDGU5 - kt * AUNDGU6 -((3*D)/(ρ*r*T)) * AUNDGU6 *(SGU6 - (ADISGU6/VGU))
        ADISGU6' = kt * ADISGU5 - kt * ADISGU6 +((3*D)/(ρ*r*T)) * AUNDGU6 *(SGU6 - (ADISGU6/VGU)) - kilGU6*ADISGU6 - kaGU*ADISGU6
        ADEGGU6' = kt*ADEGGU5 - kt*ADEGGU6 + kilGU6 * ADISGU6
        AABSGU6' = kaGU * ADISGU6

        AUNDGU7' = kt * AUNDGU6 - kt * AUNDGU7 -((3*D)/(ρ*r*T)) * AUNDGU7 *(SGU7 - (ADISGU7/VGU))
        ADISGU7' = kt * ADISGU6 - kt * ADISGU7 +((3*D)/(ρ*r*T)) * AUNDGU7 *(SGU7 - (ADISGU7/VGU)) - kilGU7*ADISGU7 - kaGU*ADISGU7
        ADEGGU7' = kt*ADEGGU6 - kt*ADEGGU7 + kilGU7 * ADISGU7
        AABSGU7' = kaGU * ADISGU7
        
        # Colon compartment
        AUNDCO' = kt * AUNDGU7 - ktCO * AUNDCO - ((3*D)/(ρ*r*T)) * AUNDCO * (SCO - (ADISCO/VCO))
        ADISCO' = kt * ADISGU7 - kt * ADISCO + ((3*D)/(ρ*r*T)) * AUNDCO * (SCO - (ADISCO/VCO)) - kilCO*ADISCO -  kaCO*ADISCO + (CP*CLI*VLI*kbil)/( Kp)
        ADEGCO' = kt * ADEGGU7 - ktCO * ADEGGU7 + kilCO * ADISCO
        AABSCO' = kaCO * ADISGU7

        #Total intestinal absorption (IA)
        AIA' = kaGU*ADISGU1 + kaGU*ADISGU2 + kaGU*ADISGU3 + kaGU*ADISGU4 + kaGU*ADISGU5 + kaGU*ADISGU6 + kaGU*ADISGU7
        
        #Somatic Compartments
        # Lungs
        CLU' = (QLU/VLU) *(CVE - (CLU*R)/(Kp))

        #Arterial blood (AR)
        CAR' = (1/VAR) * (QLU*(((CLU*R)/(Kp)) -CAR) + AIR)

        # Venous blood (VE)
        CVE' = (1/VVE) * (((QBR *CBR*R)/( Kp))  + ((QLI *CLI*R)/( Kp)) + ((QKI *CKI*R)/( Kp)) + ((QHR *CHR*R)/( Kp)) + ((QMU *CMU*R)/( Kp)) + ((QAD *CAD*R)/( Kp)) + ((QSK *CSK*R)/( Kp)) + ((QBO *CBO*R)/( Kp)) +((QTH *CTH*R)/( Kp)) - QLU * CVE + VIR)

        #Brain
        CBR' = (QBR/VBR)* (CAR - (CBR*R)/( Kp))
        #Heart
        CHR' = (QHR/VHR)* (CAR - (CHR*R)/( Kp))
        #Muscle
        CMU' = (QMU/VMU)* (CAR - (CMU*R)/( Kp))
        #Adipose
        CAD' = (QAD/VAD)* (CAR - (CAD*R)/( Kp))
        #Skin
        CSK' = (QSK/VSK)* (CAR - (CSK*R)/( Kp))
        #Bone
        CBO' = (QBO/VBO)* (CAR - (CBO*R)/( Kp))
        #Thymus
        CTH' = (QTH/VTH)* (CAR - (CTH*R)/( Kp))
        #Pancreas
        CPA' = (QPA/VPA)* (CAR - (CPA*R)/( Kp))
        #Spleen
        CSP' = (QSP/VSP)* (CAR - (CSP*R)/( Kp))

        #Stomach
        CST' = (1/VST2) * (QST*(CAR - ((CST*R)/(Kp))) + AABSST)
        #Gut
        CGU' = (1/VGU) * (QGU*(CAR - ((CGU*R)/(Kp))) + AIA)
        #Kidney
        CKI' = ((1/VKI) * QKI*(CAR - ((CKI*R)/(Kp)))) - ((CKI*Ker)/Kp)
        #Liver
        CLI' = (1/VLI)*(QHA*CAR + ((QGU *CGU*R)/( Kp)) + ((QPA *CPA*R)/( Kp)) + ((QSP *CSP*R)/( Kp)) + ((QST *CST*R)/( Kp)) - ((QLI *CLI*R)/( Kp)) - (CLI*CLint)/Kp )
    end
end