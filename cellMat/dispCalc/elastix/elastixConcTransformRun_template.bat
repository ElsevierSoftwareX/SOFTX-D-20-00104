@ECHO OFF


SET elastixPath=elastixLibPath
SET imPath=imReadPath
SET tformPath=tformReadPath
SET storePath=savePath
SET dFieldPath=%storePath%\dispField
SET regImPath=%storePath%\regIm
SET moving=%imPath%\movingImName
SET tform=%tformPath%\tformFileName



ECHO Generating the Deformation Field:
RMDIR "%dFieldPath%"
MKDIR "%dFieldPath%"
"%elastixPath%"\transformix  -def all -out "%dFieldPath%" -tp "%tform%"

ECHO Generating the Registered Image:
RMDIR "%regImPath%"
MKDIR "%regImPath%" 
"%elastixPath%"\transformix -in "%moving%" -out "%regImPath%" -tp "%tform%"s
