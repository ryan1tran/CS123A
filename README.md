# MyMSA                                     
MyMSA is a multiple sequence alignment algorithm that tries its best to compete   
against its professional counterparts (with varying levels of success)!   
   
Built entirely in C++, MyMSA is implemented using Smith-Waterman, a dynamic   
programming local alignment algorithm that is similar to Needleman-Wunsch.   

[Here is the link to the paper that utilizes this program.](https://docs.google.com/document/d/1TNdaqkgS9CMkeoG80-Iuz8-3moFuXdYUV8RzHfLvTpE/edit?usp=sharing "Bioinformatics Paper")

## Information
This project contains two directories:   
	- source   
	- application   

Within the 'source' directory is all the code used in this project. If you wish to   
compile it yourself, you may do so using either an IDE or through the command line.   

However, to save you time, the compiled '.exe' product is in the 'application' folder,   
along with the FASTA data being used for this project.   
   
   
## Running the Application
To run the application, first copy the 'application' folder somewhere onto your device.   
Once done, open up a command line (CMD for Windows and Terminal for Unix systems) and   
change directories to the application folder. Once inside, you will start the application   
by running the appropriate command:   
	- For Windows:	`.\MyMSA.exe data.fasta`   
	- For Unix:	`./MyMSA.exe data.fasta`   
   
If you would like to use your own FASTA data, replace 'data.fasta' with the other file's name.   
   
Once the application is running and you see a logo appear on the command line window (as
seen below), you are all set. From there, read and follow the instructions provided to you.   
   
![Main View](main-view.png?raw=true "Title")   
