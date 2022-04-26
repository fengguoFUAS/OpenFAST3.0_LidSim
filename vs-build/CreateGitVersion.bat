@ECHO off
SET IncludeFile=..\gitVersionInfo.h

<NUL SET /p IncludeTxt=#define GIT_VERSION_INFO '> %IncludeFile%
FOR /f %%a IN ('git describe --abbrev^=8 --always --tags --dirty') DO <NUL SET /p IncludeTxt=%%a>> %IncludeFile%
git describe --abbrev^=8 --always --tags --dirty > NUL
IF %ERRORLEVEL%==0 ( ECHO '>> %IncludeFile% ) else ( ECHO Unversioned from 42a5a8196529ae0349eda6d797a79461c2c03ff0 '>> %IncludeFile% )

EXIT /B 0