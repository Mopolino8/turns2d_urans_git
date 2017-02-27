#/bin/bash

fixsizes()
{
  awk '/Hint/{a[$3]=$NF} !/Hint/{for (i in a) gsub(i,a[i])}1' $1 >_$1
  mv _$1 $1
}

remove_cb()
{
  awk '{gsub(/_CB/,""); print}' $1 > _$1
  mv _$1 $1
}

remove_cd()
{
  awk '{gsub(/_CD/,""); print}' $1 > _$1
  mv _$1 $1
}

remove_params_global_bq()
{
  echo $1  
  awk '{gsub(/PARAMS_GLOBAL_BO/,"PARAMS_GLOBAL"); print}' $1 > _$1
  mv _$1 $1
}

add_params_sensitivity()
{
  awk '/PARAMS_GLOBAL/ {print; print "      USE PARAMS_SENSITIVITY"; next} !/PARAMS_SENSITIVITY/{print}' $1 > _$1
  mv _$1 $1
}




for i in `ls *_bo.f`; do 
  fixsizes $i
  remove_cb $i
  remove_cd $i
  remove_params_global_bq $i

  echo $i
done
for i in `ls *_bo.f`; do
  add_params_sensitivity $i
done


rm -f *_cb.f *~ *_cd.f params_global_bo.f



