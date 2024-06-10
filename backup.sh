#!/bin/bash
# Note: The MAS database is stored in a Docker volume and will persist when rebuilding images.
# This means the MYSQL_ROOT_PASSWORD configuration will not be updated unless you delete the volume containing the database ('docker volume rm mas_db-data' THIS WILL DELETE ALL UPLOADED DATA). 
# If you wish to change the MySQL password after building without deleting the volume, you can rebuild you images with the new MYSQL_ROOT_PASSWORD configuration set then manually change the password by entering the sql server container ('docker exec -it mas-sql-server') and changing it through the terminal.

docker exec mas-sql-server sh -c 'exec mysqldump -uroot -p"changeme" mas' > ./mas-db-backup.sql
docker exec -w /home/daemon/MAS -u daemon mas sh -c 'tar -czkO media' > ./mas-media-backup.tar.gz

if [[ -s ./mas-db-backup.sql ]]; then
  rclone sync -P --create-empty-src-dirs ./mas-db-backup.sql gdrive:
  if [[ $? -eq 0 ]]; then
    rm -f ./mas-db-backup.sql
  fi
fi

if [[ -s ./mas-media-backup.tar.gz ]]; then
  rclone sync -P --create-empty-src-dirs ./mas-media-backup.tar.gz gdrive:
  if [[ $? -eq 0 ]]; then
    # If rclone sync was successful, remove the file
    rm -f ./mas-media-backup.tar.gz
  fi
fi

# docker exec -i mas-sql-server sh -c 'exec mysql -uroot -p"$MYSQL_ROOT_PASSWORD" mas' < /path/to/your/backup_dir/mas-db-backup.sql
# docker exec -i -u daemon mas /bin/bash -c 'rm -rf /home/daemon/MAS/media/* ; tar -C /home/daemon/MAS/media -xz <&0' < /path/to/your/backup_dir/mas-media-backup.tar.gz