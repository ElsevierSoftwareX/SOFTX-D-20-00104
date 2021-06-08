@ECHO OFF


SET elastixPath=elastixLibPath
SET imPath=imReadPath
SET paramPath=parametersReadPath
SET storePath=savePath
SET dFieldPath=%storePath%\dispField
SET regImPath=%storePath%\regIm
SET iniTransfPath=iniTransfReadPath

SET fixedMask=%imPath%\fixedMaskName
SET movingMask=%imPath%\movingMaskName
SET fixed=%imPath%\fixedImName
SET moving=%imPath%\movingImName
SET param=%paramPath%\parameterFileName
SET iniTransf=%iniTransfPath%\iniTransfFileName



ECHO Parameters: %param%
ECHO Fixed Image: %fixed%
ECHO Moving Image: %moving%
ECHO Fixed Mask: %fixedMask%
ECHO Moving Mask: %movingMask%



ECHO elastix -f %fixed% -m %moving% -fMask %fixedMask% -mMask %movingMask% -p %param% -t0 %iniTransf% -out %storePath% 
"%elastixPath%"\elastix -f "%fixed%" -m "%moving%" -fMask "%fixedMask%" -mMask "%movingMask%" -p "%param%" -t0 "%iniTransf%" -out "%storePath%" 


ECHO Generating the Deformation Field:
RMDIR "%dFieldPath%"
MKDIR "%dFieldPath%"
"%elastixPath%"\transformix  -def all -out "%dFieldPath%" -tp "%storePath%\TransformParameters.0.txt"

ECHO Generating the Registered Image:
RMDIR "%regImPath%"
MKDIR "%regImPath%"
"%elastixPath%"\transformix -in "%moving%" -out "%regImPath%" -tp "%storePath%\TransformParameters.0.txt"


