--fill in map-verison 1 data
insert into "DTREADOUTMAPPING"("IOV_VALUE_ID", "CELL_MAP_VERSION", "ROB_MAP_VERSION") values(0,'cmssw_ROB','cmssw_ROS');
insert into "DTREADOUTCONNECTION" ("CONNECTION_ID", "IOV_VALUE_ID", "CELL", "CHANNEL", "DDU", "LAYER", "ROB", "ROS","SECTOR", "SUPERLAYER", "STATION", "TDC", "WHEEL" ) values(0,0,11,21,31,41,1,2,71,81,91,10,11);
insert into "DTREADOUTCONNECTION" ("CONNECTION_ID", "IOV_VALUE_ID", "CELL", "CHANNEL", "DDU", "LAYER", "ROB", "ROS","SECTOR", "SUPERLAYER", "STATION", "TDC", "WHEEL" ) values(1,0,11,21,31,41,1,1,71,81,91,10,11);
--fill in map-version 2 data
insert into "DTREADOUTMAPPING"("IOV_VALUE_ID", "CELL_MAP_VERSION", "ROB_MAP_VERSION") values(1,'my_ROB','my_ROS');
insert into "DTREADOUTCONNECTION" ("CONNECTION_ID", "IOV_VALUE_ID", "CELL", "CHANNEL", "DDU", "LAYER", "ROB", "ROS","SECTOR", "SUPERLAYER", "STATION", "TDC", "WHEEL" ) values(0,1,11,21,31,41,2,2,71,81,91,17,19);
insert into "DTREADOUTCONNECTION" ("CONNECTION_ID", "IOV_VALUE_ID", "CELL", "CHANNEL", "DDU", "LAYER", "ROB", "ROS","SECTOR", "SUPERLAYER", "STATION", "TDC", "WHEEL" ) values(1,1,11,21,31,41,2,1,71,81,91,17,19);