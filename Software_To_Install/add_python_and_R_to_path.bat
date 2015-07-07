REM This command adds python and R to the PATH. Assumes R-3.2.0 is installed based on the version distributed with SPARTA for Windows
setx path "C:\python27;C:\python27\Scripts\;C:\Program Files\R\R-3.2.0\bin\i386;C:\Program Files (x86)\GnuWin32\bin"
REM Set the environment variable R_LIBS_USER permanently to point to local R library to install edgeR
IF EXIST %userprofile%\Desktop\SPARTA_Windows-master (setx R_LIBS_USER %userprofile%\Desktop\SPARTA_Windows-master\R_local\library_local) ELSE (setx R_LIBS_USER %userprofile%\Desktop\SPARTA_Windows\R_local\library_local)