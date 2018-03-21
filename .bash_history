exit
ll
exit
ll
tar -xvf GenomeAnalysisTK-3.3-0.tar.bz2 
java -jar GenomeAnalysisTK.jar 
java -jar GenomeAnalysisTK.jar -T
java -jar GenomeAnalysisTK.jar -h
java -jar GenomeAnalysisTK.jar -T
java -jar GenomeAnalysisTK.jar -h
exit
mount -t nfs 130.238.36.254:/backup/ /mnt/nfs/
ping 130.238.36.254
ifconfig 
exit
make install
rm -rf /usr/local/bin/vcftools/
exit
ll
rm -rf /usr/local/bin/vcftools 
rm -rf /usr/local/bin/vcf*
git clone https://github.com/vcftools/vcftools.git
cd vcftools
./autogen.sh
./configure
make install
vcftools 
vcftools --vcf 1001genomes_snp-short-indel_only_ACGTN.vcf --min-alleles 2 --max-alleles 2 --maf 0.03 --out filter.d --recode --remove-indels
cd ../
rm -rf vcftools --vcf 1001genomes_snp-short-indel_only_ACGTN.vcf --min-alleles 2 --max-alleles 2 --maf 0.03 --out filter.d --recode --remove-in
rm -rf vcftools/
vcftools --vcf 1001genomes_snp-short-indel_only_ACGTN.vcf --min-alleles 2 --max-alleles 2 --maf 0.03 --out filter.d --recode --remove-indels
vcftools --gzvcf ./1001genomes_snp-short-indel_only_ACGTN.vcf --min-alleles 2 --max-alleles 2 --maf 0.03 --out filter.d --recode --remove-indels
vcftools --gzvcf ./1001genomes_snp-short-indel_only_ACGTN.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.03 --out filter.d --recode --remove-indels
java -jar ~/soft/beagle.27Jun16.b16.jar gt=./filter.d.recode.vcf out=./impute1001g_higher-160914.vcf nthreads=20
exit
mount 130.238.46.233:/volume1/dip /mnt
ping 130.238.46.233
ping www.google.se
ping 130.238.46.233
ssh yanjun@130.238.46.233
hostname
ping 130.238.46.233
mount 130.238.46.233:/volume1/dip /mnt
ping 130.238.46.233
yum-config-manager --enable
yum-config-manager --enable /home/yanjun/soft/
exit
man scp
scp -r yanjun@130.238.46.233:/volume1/home/dip.backup
scp -r /home/yanjun/test.txt yanjun@130.238.46.233:/volume1/home/dip.backup
scp -r /home/yanjun/test.txt yanjun@130.238.46.233:/volume1/home/
ssh yanjun@130.238.46.233
scp -r /home/yanjun/test.txt yanjun@130.238.46.233:/var/services/homes/yanjun
scp -rv /home/yanjun/test.txt yanjun@130.238.46.233:/var/services/homes/yanjun
scp -r /home/zheya/* yanjun@130.238.46.233:/var/services/homes/yanjun
ll /home/zheya/
cd /home/
sudo rm -rf zheya/
exit
mount.davfs http://130.238.46.233:5005 /mnt/
mount.davfs -h
mount -t davfs http://130.238.46.233:5005 /mnt/ uid=yanjun
mount.davfs http://130.238.46.233:5005 /mnt/ uid=yanjun
mount.davfs http://130.238.46.233:5005 /mnt/ 
vim /etc/davfs2/davfs2.conf 
mount.davfs http://130.238.46.233:5005 /mnt/ 
df -h
unmount http://130.238.46.233:5005
umount http://130.238.46.233:5005
apt-get install neon -y 
apt-get install "neon"
lsb_release -a

mount.davfs http://130.238.46.233:5005 /mnt/ 
l /mnt/lost+found/
ls /mnt/
umount /mnt/
apt-get install cadaver
mount.davfs http://130.238.46.233:5005 /mnt/ 
mount -t davfs http://130.238.46.233 /mnt/
cadaver http://130.238.46.233
cadaver http://130.238.46.233:5005
mkdir /usr/share/ca-certificates/extra
gzip /home/yanjun/archive\ \(1\).zip 
ll /home/yanjun/
ls /home/yanjun/
 ping http://130.238.46.233:5000
 ping http://130.238.46.233:5005
 ping http://130.238.46.233:5002
 ping http://130.238.46.233:5001
apt-get instll ipkg
apt-get install ipkg
cd /tmp
mkdir davfs
cd davfs
wget http://download.savannah.gnu.org/releases/davfs2/davfs2-1.4.6.tar.gz
tar xzvf davfs2-1.4.6.tar.gz
cd davfs2-1.4.6
./configure --prefix=/opt --with-neon=/opt --with-ssl=openssl
cd /tmp/
wget http://www.webdav.org/neon/neon-0.25.4.tar.gz
tar xvzf neon-0.25.4.tar.gz
cd neon-0.25.4
LIBS="-ldl" ./configure --prefix=/usr/local/neon --with-ssl=openssl --with-libxml2 --with-libs=/usr/local/openssl/:/usr/local/libxml2/ --without-zlib
apt-get install libssl-dev
LIBS="-ldl" ./configure --prefix=/usr/local/neon --with-ssl=openssl --with-libxml2 --with-libs=/usr/local/openssl/:/usr/local/libxml2/ --without-zlib
apt-get install libxml2
LIBS="-ldl" ./configure --prefix=/usr/local/neon --with-ssl=openssl --with-libxml2 --with-libs=/usr/local/openssl/:/usr/local/libxml2/ --without-zlib
apt-get install xml2
LIBS="-ldl" ./configure --prefix=/usr/local/neon --with-ssl=openssl --with-libxml2 --with-libs=/usr/local/openssl/:/usr/local/libxml2/ --without-zlib
apt-get install libxml2-devel
apt-get install libxml2
yum install libxml2-devel
aptitude serch libxml2
aptitude serch libxml2
aptitude search libxml2
apt-get install 
apt-get install libxml2-dev
apt-get install  libxml2-utils
LIBS="-ldl" ./configure --prefix=/usr/local/neon --with-ssl=openssl --with-libxml2 --with-libs=/usr/local/openssl/:/usr/local/libxml2/ --without-zlib
make
make install
mount.davfs http://130.238.46.233:5005 /mnt/
mount.davfs http://130.238.46.233:5005 /mnt/
cp /home/yanjun/archive\ \(1\)/* /etc/davfs2/certs/
cp /home/yanjun/archive\ \(2\)/* /etc/davfs2/certs/
mount.davfs http://130.238.46.233:5005 /mnt/
man davfs2.conf
vim /etc/davfs2/davfs2.conf 
mount.davfs http://130.238.46.233:5005 /mnt/
df -h
umount /mnt
mount.davfs http://130.238.46.233:5005 /mnt/
mount.davfs http://130.238.46.233 /mnt/
umount /mnt/
mount.davfs http://130.238.46.233 /mnt/
ls /mnt/
umount /mnt
mount.davfs http://130.238.46.233:5000 /mnt/
mount.davfs https://130.238.46.233:50001 /mnt/
mount.davfs https://130.238.46.233:5001 /mnt/
cd /etc/davfs2/certs/
ll
ls
less -S syno-ca-cert.pem
less -S syno-ca-privkey.pem 
vim /etc/davfs2/davfs2.conf 
mount.davfs https://130.238.46.233:5001 /mnt/
umount /mnt
vim /etc/davfs2/davfs2.conf 
mount.davfs https://130.238.46.233:5001 /mnt/
vim /etc/davfs2/davfs2.conf 
mount.davfs https://130.238.46.233:5001 /mnt/
vim /etc/davfs2/davfs2.conf 
openssl s_client -connect https://130.238.46.233:5001 -showcerts
openssl s_client -connect https://130.238.46.233:5001 -showcerts </dev/null 2>/dev/null
openssl s_client -connect https://130.238.46.233:5001 -showcerts </dev/null 2>/dev/null
openssl s_client -connect https://130.238.46.233:5001 -showcerts 
openssl s_client -connect https://130.238.46.233:5001 -showcerts  </dev/null 2>/dev/null | openssl x509 -outform PEM > certificate.pem
ll
ls
less -s certificate.pem 
openssl s_client -connect https://130.238.46.233:443 -showcerts  </dev/null 2>/dev/null | openssl x509 -outform PEM > certificate.pem
less -s certificate.pem 
pt install ca-certificates
apt install ca-certificates
mount.davfs https://130.238.46.233:5001 /mnt/
df -h
ls /mnt
ls -al /mnt
mount.davfs dav://130.238.46.233:5005 /mnt/
umount /mnt
mount.davfs http://130.238.46.233:5005 /mnt/
umount /mnt
ls
less certificate.pem
less cert.pem
vim /etc/davfs2/davfs2.conf 
mount.davfs http://130.238.46.233:5005 /mnt/
mount.davfs http://130.238.46.233:5000 /mnt/
ll
ls
rm *
rm -rf *
ls
cp /home/yanjun/archive\ \(4\)/* ./
ll
ls
vim /etc/davfs2/davfs2.conf 
mount.davfs http://130.238.46.233:5000 /mnt/
df -h
umount /mnt
mount.davfs http://130.238.46.233:5000 /mnt/
df -h
df -h
mount.davfs https://130.238.46.233:5001 /mnt/
l
ls
vim /etc/davfs2/davfs2.conf 
mount.davfs https://130.238.46.233:5001 /mnt/
ifconfig
mount 130.238.46.233:/volume1/BBG /mnt/
mount 130.238.46.233:/volume1/BBG /mnt/
mount_nfs 130.238.46.233:/volume1/BBG /mnt/
sudo apt-get install nfs-kernel-server
sudo apt-get install nfs-common
mount
man mount
mount -t nfs 130.238.46.233:/volume1/BBG/ /mnt
mount 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs -o resvport 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs 130.238.46.233:/volume1/BBG/ /mnt
/etc/init.d/rpcbind start
mount -t nfs 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs -o resvport 130.238.46.233:/volume1/BBG/ /mnt
mount 
man mount
mount -t nfs -v 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs -o resvport -v 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs -v o nfsvers=2 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs -v  nfsvers=2 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs -v o nfsvers=2 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs -v  nfsvers=2 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs -o resvport -nfsvers=2 -v 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs -o resvport -nfsvers 2 -v 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs  -v 130.238.46.233:/volume1/BBG/ /mnt
umount /mnt
umount /mnt/
df -h
mount -t nfs -o nfsvers=2  -v 130.238.46.233:/volume1/BBG/ /mnt
mount -t nfs -o nfsvers=3  -v 130.238.46.233:/volume1/BBG/ /mnt
rpcinfo -p 
mount.nfs  -v 130.238.46.233:/volume1/BBG/ /mnt
mount.nfs -o resvport -v 130.238.46.233:/volume1/BBG/ /mnt
nslookup 192.168.1.5
nslookup 130.238.46.233
 sudo apt-get install libpam-krb5
nslookup diprotodon.imbim.uu.se
nslookup 130.238.46.233
cat /etc/fstab
man mount
mount.davfs 130.238.46.233:5005 /mnt
ping 130.238.46.233
lsb_release 
lsb_release -v
lsb_release -version
man lsb_release 
lsb_release -v
lsb_release -version
lsb_release -
lsb_release -v
dpkg -l | grep nfs-common
lsb_release -a
sudo apt-get update 4.5.0-040500rc7-generic
apport information
mount -t davfs https://130.238.46.233:5001 /mnt/
vim /etc/davfs2/davfs2.conf 
mount -t davfs -v https://130.238.46.233:5001 /mnt/
sudo a2enmod ssl
  service apache2 restart
a2ensite default-ssl
  service apache2 reload
sudo apt-get install ssl-cert
make-ssl-cert generate-default-snakeoil --force-overwrite
sudo /etc/init.d/apache2 restart
mount -t davfs -v https://130.238.46.233:5001 /mnt/
umount /mnt
mount -t davfs -v https://130.238.46.233:5001 /mnt/
apt-get install rpcbind nfs-common
rpcbind : ALL
rcpbind 130.238.46.233
apt-get install rpcbind
rcpbind 130.238.46.233
rpcbind 130.238.46.233
rpcbind 130.238.46.233
rpcbind :ALL
mount 130.238.46.233:/volume1/BBG /mnt/
mount -v 130.238.46.233:/volume1/BBG /mnt/
cat /etc/host
cat /etc/hosts
vim /etc/hosts
cat /etc/hosts
cat /etc/hosts
vim /etc/hosts
cat /etc/hosts
mount -v 130.238.46.233:/volume1/BBG /mnt/
mount.davfs -v http://130.238.46.233:5005 /mnt/
umount mnt
umount /mnt
mount.davfs -v http://130.238.46.233:5005 /mnt/
umount /mnt
vim /etc/davfs2/davfs2.conf 
vim ~/.davfs2/secrets
sudo vim ~/.davfs2/secrets
vim /etc/davfs2/davfs2.conf 
mount -t davfs http://130.238.46.233:5005 /mnt
mount -t davfs -v http://130.238.46.233:5005 /mnt
mount -t nfs -o vers=2 130.238.46.233:/volume1/BBG /mnt
mount -t nfs -v -o vers=2 130.238.46.233:/volume1/BBG /mnt
mount -t nfs -v -o vers=4 130.238.46.233:/volume1/BBG /mnt
sudo apt-get install cifs-utils
mount -t nfs -v -o vers=4 130.238.46.233:/volume1/BBG /mnt
mount.nfs -v 130.238.46.233:/volume1/BBG /mnt
exit
mount.davfs http://130.238.46.233:5005 /mnt
umount /mnt
ping wombat3.synology.me
exit
pwd
cd /tmp/
wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2
tar -xvzf fastx_toolkit-0.0.14.tar.bz2 
gzip -d fastx_toolkit-0.0.14.tar.bz2 
 $ tar -xjf libgtextutils-0.6.tar.bz2
 tar -xjf libgtextutils-0.6.tar.bz2
wget https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz
tar -xjf libgtextutils-0.7.tar.gz 
ls
gzip libgtextutils-0.7.tar.gz 
ll
ls 
gzip libgtextutils-0.7.tar.gz
gzip -d libgtextutils-0.7.tar.gz
cd libgtextutils-0.7.tar 

cd libgtextutils-0.7.tar 
ls
tar -zxvf libgtextutils-0.7.tar
wget https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz
tar -zxvf libgtextutils-0.7.ta
tar -zxvf libgtextutils-0.7.tar.gz 
cd libgtextutils-0.7
./configure
make 
make install
cd ..
ls
tar -xjf fastx_toolkit-0.0.12.tar.bz2 
tar -xjf fastx_toolkit-0.0.14.tar.bz2 
cd fastx_toolkit-0.0.14
./configure
make 
make install
exit
apt-get install nfs-common 
mount -t nfs -o proto=tcp,port=2049 130.238.46.94:/mnt/ /mnt
sudo apt-get install rpcbind nfs-common
mount -t nfs -o proto=tcp,port=2049 130.238.46.94:/mnt/ /mnt
ping 130.238.46.94
exit
exit
