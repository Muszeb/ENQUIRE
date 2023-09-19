WATCHED_PID=$({ ./code/ENQUIRE.sh -f input/textmining_config.txt >ENQUIRE.stdout 2>ENQUIRE.stderr & } && echo $!);
while [[ $(ps -p $WATCHED_PID | wc -l) -gt 1 ]]; do
	ps -T -o gid,pgid,ppid,pid,lwp,time,etime,bsdtime,%cpu,%mem,rss,comm,args --no-headers --cumulative S | grep -w $WATCHED_PID\| | sed 's/\s\+/ \"/12' | sed 's/$/\"/'
	#ps -L -F -l -c -g $(ps -o sid= -p $WATCHED_PID) --headers | grep $WATCHED_PID #| sed 's/\s\+/ \"/10' | sed 's/$/\"/'
	sleep 1 
done
#
# header
#head="pgid\tppid\tpid\ttime\tetime\tbsdtime\t%cpu\t%mem\trss\tcomm\targs\n"
#echo -e $head
#sed -i "1s/^/${head}/" ENQUIRE_stats2.txt 