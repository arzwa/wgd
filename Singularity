Bootstrap: docker
From: ubuntu

%runscript
	wgd --help

%help
	wgd --help

%labels
	AUTHOR arzwa@psb.vib-ugent.be

%environment
	# click was complaining about this
	export LC_ALL=C.UTF-8
	export LANG=C.UTF-8

%post
	# install python, git, etc.
	apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -yq install python3-pip python3-tk git wget \
	    build-essential mcl ncbi-blast+ muscle mafft prank fasttree phyml paml

	# set an alias for fasttree
	ln -s /usr/bin/fasttree /usr/bin/FastTree

	# get wgd
	git clone https://github.com/arzwa/wgd.git
	cd wgd

	# install wgd
	pip3 install .
