@ECHO OFF


SET elastixPath=elastixLibPath
SET imPath=imReadPath
SET paramPath=parametersReadPath
SET storePath=savePath
SET fixed=%imPath%\fixedImName
SET moving=%imPath%\movingImName
SET param=%paramPath%\parameterFileName


ECHO Fixed Image: %fixed%
ECHO Moving Image: %moving%
ECHO Parameters: %param%


ECHO elastix -f %fixed% -m %moving% -p %param% -out %storePath% 
"%elastixPath%"\elastix -f "%fixed%" -m "%moving%" -p "%param%" -out "%storePath%"