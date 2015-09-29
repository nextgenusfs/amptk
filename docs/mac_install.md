####Mac OSX Install Instructions####
___

Note: there are many ways to install these dependencies, this is one way.

#####1) Need to install Xcode (Developer Tools that Apple doesn’t install natively but supports):#####
```xcode-select --install```

#####2) Install HomeBrew (copy and paste this into terminal):#####
```ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)”```

Then setup homebrew: type `brew doctor`, then type: `brew tap homebrew/science`

#####3) Install some tools

```
#install python via homebrew
brew install python

#use pip to install biopython and natsort, might require sudo
pip install biopython natsort  
```

#####4) Now download this repository:

`git clone https://github.com/nextgenusfs/ufits`

And then you will need to add the following to your `~/.bash_profile` or always include the entire path to the scripts at runtime.

`export PATH="/location/of/packages/ufits:$PATH"`

#####5) You will also need to install USEARCH8 - get it here: http://www.drive5.com/usearch/download.html.  One way to make the program executable and move into your path:

```
#make executable
sudo chmod +x /path/to/usearch8.0.1623_i86osx32
```

```
#create softlink to folder in $PATH
sudo ln -s /path/to/usearch8.0.1623_i86osx32 /usr/local/bin/usearch8
```

