ensure docker is running  
git clone https://github.com/ekazasi/Kazasi_Eftychia_Assignment.git  
cd Kazasi_Eftychia_Assignment  
docker build -t bio_final_image .  
docker run -it bio_final_image  
./all_script.sh config.yaml
