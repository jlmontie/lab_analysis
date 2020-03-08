for RUNDIR in /srv/idbydna-group3/results/idbd_rnd_v2/200209_NB551543_0214_AH7FHKBGXF/ \
/srv/idbydna-group3/results/idbd_rnd_v2/200214_NB551702_0118_AH7FC2BGXF/ \
/srv/idbydna-group3/results/idbd_rnd_v2/200215_NB551543_0219_AH7CJ5BGXF/ \
/srv/idbydna-group3/results/idbd_rnd_v2/200218_NB551543_0221_AH7FHJBGXF/

do
    python ./../get_nr_from_rundir.py $RUNDIR -g 1578
done