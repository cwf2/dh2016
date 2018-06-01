# Install prerequisites on Vagrant VM
# using ubuntu/trusty64

sudo -u vagrant cp /vagrant/vagrant/screenrc /home/vagrant/.screenrc
sudo -u vagrant cp /vagrant/vagrant/vimrc /home/vagrant/.vimrc

apt-get install -y \
  git \
  htop \
  libxml2-dev \
  libgsl0-dev

# install R

apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list

apt-get update
apt-get install -y r-base-dev

# install R packages

Rscript --vanilla /vagrant/vagrant/install.R

