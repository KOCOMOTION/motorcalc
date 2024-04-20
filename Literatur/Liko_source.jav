
U = 0.0;                            // Voltage
I = 0.0;                            // Current
Rm = 0.0;                           // Resistance
Ui = 0.0;                           // induced Voltage
nM = 0.0;                           // speed?
nL = 0.0;                           // load-speed?
M = 0.0;                            // torque
MM = 0.0;                           // max. torque?
ML = 0.0;                           // load torque?
Io = 0.0;                           // no-load current
k = 0.0;                            // motor constant 
backEMF = 0.0;                      // backEMF constant
ratio = 0.0;                        // gear ratio   
gear_Efficiency = 0.0;              // gear efficiency    
system_Efficiency = 0.0;
Mf = 0.0;                           
PL = 0.0;
PM = 0.0;
T = 0.0;
regulation = 0.0;
nMCharacteristic = 0.0;
maxMM = null;
PME = 0.0;
T_ambient = 0.0;
R = 0.0;
maxI = null;
maxT = null;
Rth1 = null;
Rth2 = null;
nGMax = null;
MGMax = null;
stallTorque = 0.0;
stallCurrent = 0.0;
n_o = 0.0;


maxI = getMotor(context).getMaxCurrent(context);
nMMax = getMotor(context).getMaxSpeed(context);
maxT = getMotor(context).getMaxCoilTemperature(context);

nGMax = getGearbox(context).getMaxSpeed(context);
MGMax = getCalculation().getRatioMaxTorque(context);

double temperaturkoeffizientKupfer = 0.004;

R = getMotor(context).getResistance(context);
k = getMotor(context).getTorqueConstant(context);
backEMF = (k * 1000 * 2 * Math.PI) / 60;
Io = getMotor(context).getNoLoadCurrent(context);
T_ambient = getCalculation().getAmbientTemp(context);
ratio = getCalculation().getRatio(context);
gear_Efficiency = getCalculation().getRatioEfficiency(context);


Mf = Io * k;

switch (getCalculation().getLastIndexOfTorqueEdit(context)) {
    case Calculation.LAST_IndexOfTorqueEdit_MotorCurrent:
        I = getCalculation().getLastValueOfTorqueEdit(context);
        M = I * k;
        MM = M - Mf;
        ML = MM * ratio * gear_Efficiency;
        break;
    case Calculation.LAST_IndexOfTorqueEdit_MotorTorque:
        MM = getCalculation().getLastValueOfTorqueEdit(context);
        M = MM + Mf;
        I = M / k;
        ML = MM * ratio * gear_Efficiency;
        break;
    case Calculation.LAST_IndexOfTorqueEdit_LoadTorque:
        ML = getCalculation().getLastValueOfTorqueEdit(context);
        MM = ML / (ratio * gear_Efficiency);
        M = MM + Mf;
        I = M / k;
        break;
}

double T22 = 22; // temperature of Resistance
T =  T_ambient;
if( getMotor( context ).getMotorTechnology( context ).getAnkerType( context ) == MotorTechnology.AnkerType_Coreless ){
    double Rthg = 0.0;

    dutiCycleTime = getCalculation().getDutiCycleTime(context);
    Rth1 = getMotor(context).getThermalResistanceRotorBody(context);
    tth1 = getMotor(context).getThermalTimeConstantRotor(context);
    Rth2 = getMotor(context).getThermalResistanceBodyAmbient(context);
    tth2 = getMotor(context).getThermalTimeConstantBody(context);

    if ((dutiCycleTime != null) && (tth1 != null)) {
        Rthg += ( Rth1 * ( 1 - Math.exp((-1) * dutiCycleTime / tth1 )));
    } else {
        Rthg += Rth1;
    }

    if ((dutiCycleTime != null) && (tth2 != null)) {
        Rthg += ( Rth2 * ( 1 - Math.exp((-1) * dutiCycleTime / tth2 )));
    } else {
        Rthg += Rth2;
    }

    double a = 1.0 - Rthg * R * I * I * temperaturkoeffizientKupfer;
    
    if (a <= 0.0) {
        a = Math.sqrt( 1.0 / (Rthg * R * temperaturkoeffizientKupfer));
        try {
            format = new LDecimalFormat();
            format.setMinScale(0);
            format.setMaxScale(4);
            format.setRoundAfterDigits(LikoGlobalDelegation.DEFAULT_RoundAfterDigits);
            str = format.valueToString( a );
        } catch (LFormatException exception) {
            str = "";
        }

        list.add(LikoGlobalDelegation.getInstance().getResourceBundle().getString("CalculationController.startCalculation.TemperatureInfinite1"));
        list.add(
            String.format(
            LikoGlobalDelegation.getInstance().getResourceBundle().getString("CalculationController.startCalculation.TemperatureInfinite2"),
            str));
        T = 999999999;
    } else {
        T = (Rthg * R * I * I * (1.0 - temperaturkoeffizientKupfer * T22) + T_ambient) / a;
    }
}

Rm = R * (1.0 + temperaturkoeffizientKupfer * (T - T22));
regulation = Rm / (k * k);
nMCharacteristic = ( regulation * 60 ) / ( 2 * Math.PI );

double URm = I * Rm;

switch (getCalculation().getLastIndexOfSpeedEdit(context)) {
    case Calculation.LAST_IndexOfSpeedEdit_MotorVoltage:
        U = getCalculation().getLastValueOfSpeedEdit(context);
        Ui = U - URm;
        nM = (Ui * 1000) / backEMF;
        nL = nM / ratio;
        break;
    case Calculation.LAST_IndexOfSpeedEdit_MotorSpeed:
        nM = getCalculation().getLastValueOfSpeedEdit(context);
        Ui = (nM * backEMF) / 1000;
        U = Ui + URm;
        nL = nM / ratio;
        break;
    case Calculation.LAST_IndexOfSpeedEdit_LoadSpeed:
        nL = getCalculation().getLastValueOfSpeedEdit(context);
        nM = nL * ratio;
        Ui = (nM * backEMF) / 1000;
        U = Ui + URm;
        break;
}

if( getMotor( context ).getMotorTechnology( context ).getAnkerType( context ) == MotorTechnology.AnkerType_Core ){
    currentAtMaxEfficiency = Math.sqrt( Io  * (( U / Rm )  - Io ));
}

stallTorque = ( U / Rm ) * k - Mf;
stallCurrent = ( U / R );
n_o = (( U - ( Io * R )) * 60 ) / ( 2 * Math.PI * k );

double wM = nM * 2 * Math.PI / 60;
PM = MM * wM;
PME = U * I;
PL = ML * wM / ratio;
system_Efficiency = PL / PME;

maxI = getMotor(context).getMaxCurrent(context);
/*
if(( maxI == null ) && ( maxT != null )){
double b = Rthg * R * ( 1 + temperaturkoeffizientKupfer * ( maxT.doubleValue() - T22 ));
maxI = new Double(
Math.sqrt(( maxT.doubleValue() - T_ambient ) / b )
);
}
    */
if (maxI != null) {
    maxMM = (maxI - Io) * k;
}

if (I < Io) {
    try {
        format = new LDecimalFormat();
        format.setMinScale(0);
        format.setMaxScale(4);
        format.setRoundAfterDigits(LikoGlobalDelegation.DEFAULT_RoundAfterDigits);
        str = format.valueToString(Io);
    } catch (LFormatException exception) {
        str = "";
    }
    list.add(
        String.format(
            LikoGlobalDelegation.getInstance().getResourceBundle().getString("CalculationController.startCalculation.CurrentLowerNoLoadCurrent"),
            str
        )
    );
    noError = false;
}

if (U <= URm) {
    try {
        format = new LDecimalFormat();
        format.setMinScale(0);
        format.setMaxScale(4);
        format.setRoundAfterDigits(LikoGlobalDelegation.DEFAULT_RoundAfterDigits);
        str = format.valueToString(URm);
    } catch (LFormatException exception) {
        str = "";
    }
    list.add(
        String.format(
            LikoGlobalDelegation.getInstance().getResourceBundle().getString("CalculationController.startCalculation.VoltageHigherStartingVoltage"),
            str
        )
    );
    noError = false;
}

if ((maxI != null) && (maxI < I)) {
    list.add(LikoGlobalDelegation.getInstance().getResourceBundle().getString("CalculationController.startCalculation.CurrentHigherMaxCurrent"));
}

if ((nMMax != null) && (nMMax < nM)) {
    list.add(LikoGlobalDelegation.getInstance().getResourceBundle().getString("CalculationController.startCalculation.MotorSpeedHigherMaxSpeed"));
}

if ((nGMax != null) && (nGMax < nM)) {
    list.add(LikoGlobalDelegation.getInstance().getResourceBundle().getString("CalculationController.startCalculation.GearBoxSpeedHigherMaxSpeed"));
}

if ((maxT != null) && (maxT < T)) {
    list.add(LikoGlobalDelegation.getInstance().getResourceBundle().getString("CalculationController.startCalculation.CoilTempHigherMaxCoilTemp"));
}

