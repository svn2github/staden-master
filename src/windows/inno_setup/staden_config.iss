; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{FFFFFFFF-FF03-4003-BFEA-0008028F89FC}
AppName=Staden Package
AppVersion=2.0.0b10
;AppVerName=Staden Package 2.0.0b10
AppPublisher=Wellcome Trust Sanger Institute
AppPublisherURL=http://staden.sourceforge.net
AppSupportURL=http://staden.sourceforge.net
AppUpdatesURL=http://staden.sourceforge.net
ChangesAssociations=yes
DefaultDirName={pf}\Staden Package
;DisableDirPage=yes
DefaultGroupName=Staden Package
DisableProgramGroupPage=yes
LicenseFile=C:\staden-build-{%OS_SIZE}\LICENCE.txt
OutputDir=C:\staden-build-{%OS_SIZE}\windows\inno_setup
OutputBaseFilename=staden_setup
Compression=lzma
;Compression=zip/1
SolidCompression=yes
; none => install as normal user if you have no admin rights, otherwise install as administrator
PrivilegesRequired=none
WizardImageFile=c:\staden-build-{%OS_SIZE}\windows\inno_setup\splash.bmp

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Files]
Source: "C:\staden-inst-{%OS_SIZE}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

;[Tasks]
;Name: associate; Description: "&Associate files"; GroupDescription: "Other tasks:"; Flags: unchecked

[Icons]
Name: "{group}\Gap5"; Filename: "{app}\bin\Gap5.exe"
Name: "{group}\Gap4"; Filename: "{app}\bin\Gap.exe"
Name: "{group}\Trev"; Filename: "{app}\bin\Trev.exe";
Name: "{group}\Spin"; Filename: "{app}\bin\Spin.exe"
Name: "{group}\Pregap4"; Filename: "{app}\bin\Pregap4.exe"
Name: "{group}\URL - Home Page"; Filename: "http://staden.sourceforge.net"
Name: "{group}\URL - Local documentation"; Filename: "{app}\share\doc\staden\index.html"
Name: "{group}\{cm:UninstallProgram,Staden Package}"; Filename: "{uninstallexe}"

[Registry]
; File types for admin
Root: HKLM; Subkey: "Software\Classes\.ztr"; ValueType: string; ValueName: ""; ValueData: "ZTRTraceFile"; Flags: uninsdeletevalue; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\ZTRTraceFile"; ValueType: string; ValueName: ""; ValueData: "ZTR Trace File"; Flags: uninsdeletekey; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\ZTRTraceFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\trev.exe,0"; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\ZTRTraceFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\trev.exe"" ""%1"""; Check: IsAdmin

Root: HKLM; Subkey: "Software\Classes\.ab1"; ValueType: string; ValueName: ""; ValueData: "AB1TraceFile"; Flags: uninsdeletevalue; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\AB1TraceFile"; ValueType: string; ValueName: ""; ValueData: "AB1 Trace File"; Flags: uninsdeletekey; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\AB1TraceFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\trev.exe,0"; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\AB1TraceFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\trev.exe"" ""%1"""; Check: IsAdmin

Root: HKLM; Subkey: "Software\Classes\.exp"; ValueType: string; ValueName: ""; ValueData: "ExperimentFile"; Flags: uninsdeletevalue; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\ExperimentFile"; ValueType: string; ValueName: ""; ValueData: "Experiment File File"; Flags: uninsdeletekey; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\ExperimentFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\trev.exe,0"; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\ExperimentFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\trev.exe"" ""%1"""; Check: IsAdmin

Root: HKLM; Subkey: "Software\Classes\.g5d"; ValueType: string; ValueName: ""; ValueData: "Gap5DBFile"; Flags: uninsdeletevalue; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Gap5DBFile"; ValueType: string; ValueName: ""; ValueData: "Gap5 Database File"; Flags: uninsdeletekey; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Gap5DBFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\gap5.exe,0"; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Gap5DBFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\gap5.exe"" ""%1"""; Check: IsAdmin

Root: HKLM; Subkey: "Software\Classes\.g5x"; ValueType: string; ValueName: ""; ValueData: "Gap5AuxFile"; Flags: uninsdeletevalue; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Gap5AuxFile"; ValueType: string; ValueName: ""; ValueData: "Gap5 Index File"; Flags: uninsdeletekey; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Gap5AuxFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\gap5.exe,0"; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Gap5AuxFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\gap5.exe"" ""%1"""; Check: IsAdmin

Root: HKLM; Subkey: "Software\Classes\.aux"; ValueType: string; ValueName: ""; ValueData: "Gap4AuxFile"; Flags: uninsdeletevalue; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Gap4AuxFile"; ValueType: string; ValueName: ""; ValueData: "Gap4 Index File"; Flags: uninsdeletekey; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Gap4AuxFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\gap.exe,0"; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Gap4AuxFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\gap.exe"" ""%1"""; Check: IsAdmin

Root: HKLM; Subkey: "Software\Classes\.pg4"; ValueType: string; ValueName: ""; ValueData: "Pregap4ConfigFile"; Flags: uninsdeletevalue; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Pregap4ConfigFile"; ValueType: string; ValueName: ""; ValueData: "Pregap4 Configuration File"; Flags: uninsdeletekey; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Pregap4ConfigFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\preggap4.exe,0"; Check: IsAdmin
Root: HKLM; Subkey: "Software\Classes\Pregap4ConfigFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\pregap4.exe"" ""%1"""; Check: IsAdmin

; File types for non-admin
Root: HKCU; Subkey: "Software\Classes\.ztr"; ValueType: string; ValueName: ""; ValueData: "ZTRTraceFile"; Flags: uninsdeletevalue; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\ZTRTraceFile"; ValueType: string; ValueName: ""; ValueData: "ZTR Trace File"; Flags: uninsdeletekey; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\ZTRTraceFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\trev.exe,0"; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\ZTRTraceFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\trev.exe"" ""%1"""; Check: IsNotAdmin

Root: HKCU; Subkey: "Software\Classes\.ab1"; ValueType: string; ValueName: ""; ValueData: "AB1TraceFile"; Flags: uninsdeletevalue; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\AB1TraceFile"; ValueType: string; ValueName: ""; ValueData: "AB1 Trace File"; Flags: uninsdeletekey; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\AB1TraceFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\trev.exe,0"; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\AB1TraceFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\trev.exe"" ""%1"""; Check: IsNotAdmin

Root: HKCU; Subkey: "Software\Classes\.exp"; ValueType: string; ValueName: ""; ValueData: "ExperimentFile"; Flags: uninsdeletevalue; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\ExperimentFile"; ValueType: string; ValueName: ""; ValueData: "Experiment File File"; Flags: uninsdeletekey; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\ExperimentFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\trev.exe,0"; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\ExperimentFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\trev.exe"" ""%1"""; Check: IsNotAdmin

Root: HKCU; Subkey: "Software\Classes\.g5d"; ValueType: string; ValueName: ""; ValueData: "Gap5DBFile"; Flags: uninsdeletevalue; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Gap5DBFile"; ValueType: string; ValueName: ""; ValueData: "Gap5 Database File"; Flags: uninsdeletekey; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Gap5DBFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\gap5.exe,0"; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Gap5DBFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\gap5.exe"" ""%1"""; Check: IsNotAdmin

Root: HKCU; Subkey: "Software\Classes\.g5x"; ValueType: string; ValueName: ""; ValueData: "Gap5AuxFile"; Flags: uninsdeletevalue; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Gap5AuxFile"; ValueType: string; ValueName: ""; ValueData: "Gap5 Index File"; Flags: uninsdeletekey; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Gap5AuxFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\gap5.exe,0"; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Gap5AuxFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\gap5.exe"" ""%1"""; Check: IsNotAdmin

Root: HKCU; Subkey: "Software\Classes\.aux"; ValueType: string; ValueName: ""; ValueData: "Gap4AuxFile"; Flags: uninsdeletevalue; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Gap4AuxFile"; ValueType: string; ValueName: ""; ValueData: "Gap4 Index File"; Flags: uninsdeletekey; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Gap4AuxFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\gap.exe,0"; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Gap4AuxFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\gap.exe"" ""%1"""; Check: IsNotAdmin

Root: HKCU; Subkey: "Software\Classes\.pg4"; ValueType: string; ValueName: ""; ValueData: "Pregap4ConfigFile"; Flags: uninsdeletevalue; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Pregap4ConfigFile"; ValueType: string; ValueName: ""; ValueData: "Pregap4 Configuration File"; Flags: uninsdeletekey; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Pregap4ConfigFile\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\bin\preggap4.exe,0"; Check: IsNotAdmin
Root: HKCU; Subkey: "Software\Classes\Pregap4ConfigFile\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\bin\pregap4.exe"" ""%1"""; Check: IsNotAdmin

[Code]
function IsAdmin(): Boolean;
begin
  Result := IsAdminLoggedOn or IsPowerUserLoggedOn;
end;

function IsNotAdmin(): Boolean;
begin
  Result := not (IsAdminLoggedOn or IsPowerUserLoggedOn);
end;