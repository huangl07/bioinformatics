# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

export package=/home/env/package
export prefix=/home/env/publib
export bin=$prefix/bin/
export biotools=/home/env/biotools/
export CFLAGS="-fPIC"
export GOROOT=$prefix/go
export PATH=$GOROOT/bin:$prefix/bin:/usr/bin/
export LD_LIBRARY_PATH=$prefix/lib/:$prefix/lib64/
export LIBRARY_PATH=$prefix/lib/:$prefix/lib64/
export C_INCLUDE_PATH=$prefix/include/
export CPLUS_INCLUDE_PATH=$prefix/include/
export PKG_CONFIG_PATH=$prefix/lib/pkgconfig:$prefix/lib64/pkgconfig/
export MANPATH=$prefix/share/man:$MANPATH
export JAVA_HOME=$prefix/java
export PATH=$JAVA_HOME/bin:$PATH 
export CLASSPATH=.:$JAVA_HOME/lib/dt.jar:$JAVA_HOME/lib/tools.jar 
R_HOME=$prefix/lib64/R
export PATH=$R_HOME/bin/:$PATH
export LD_LIBRARY_PATH=$R_HOME/lib/:$LD_LIBRARY_PATH
export ROOTSYS=$prefix/root/
export LD_LIBRARY_PATH=$ROOTSYS/lib/:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin/:$PATH
red='\e[01;31m'
cyan='\e[01;36m'
blue='\e[01;34m'
yellow='\e[01;33m'
purple='\e[01;35m'
green='\e[01;32m'
white='\e[0;37m'
PS1="$red\u $yellow: $green\h $yellow: $blue\t $yellow: $purple\d $yellow: $cyan\w $white\n$"
