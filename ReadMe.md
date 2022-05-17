The file structure of the project is complicated because it contains 4 separate tutorial projects provided by ST Micro. 
Furthermore, even though the project is a C++ project, the main file is actually a C file that is recognized as a C++ file by the IDE. 
The reason for this is unknown. The project that runs the Changepoint detector is the Datalog project. 
Within this project folder, is another STm32CubeIDE folder that contains a folder structure for supplemental files to the Datalog project. 

The User folder contains C++ files that contain code for the change point algorithm. In order to load the project into the StM32CubeIDE, the root folder can be imported into the IDE as an existing project directory. 
The Import wizard will identify all of the different projects to import, but only the Datalog project needs to be selected to import the project. 
