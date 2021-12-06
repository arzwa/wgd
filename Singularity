Bootstrap: docker
From: vibpsb/i-adhore

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
	# set TMPDIR to a path defined within the container for mafft execution
	export TMPDIR="/tmp"

%post
	# install PAML from source
	apt-get install -y wget && wget http://abacus.gene.ucl.ac.uk/software/paml4.9j.tgz && \
	tar -xzf paml4.9j.tgz && cd paml4.9j/src && make -f Makefile && mv codeml /bin && cd /

	# install python, git, etc.
	apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -yq install python3-pip python3-tk git \
	    build-essential diamond-aligner mcl ncbi-blast+ muscle mafft prank fasttree phyml

	# get wgd
	git clone https://github.com/arzwa/wgd.git
	cd wgd

	# install wgd
	pip3 install .