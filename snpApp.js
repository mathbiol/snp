console.log('snpApp :-)')

// Assemble application under div id=snpApp

// 1. Get SNP textArea going
snpApp={

	disp:function(x){
		//console.log(x);
		var br = document.createElement('li');
		br.innerHTML=x; // note innerHTML allows HTML tagging, not as plain as textContent
		jmat.gId('disp').appendChild(br)},

	start:function(){
		// add text area to collect SNP string
		//ta = document.createElement('textarea');
		//ta.id='snpMatrix';
		//jmat.gId('snpApp').appendChild(ta);
		var ta = jmat.gId('snpCalc')
		ta.onclick=function(evt){
			jmat.gId('disp').innerHTML='hr';
			snpApp.disp('starting ...');
			snpApp.parseSnp(evt);
			snpApp.disp('SNPs: '+snpApp.dt.snps[0].length);
			var n=snpApp.dt.rows.length;
			snpApp.disp('Strains:');
			for (i=0;i<n;i++){snpApp.disp(i+'/'+n+' - '+snpApp.dt.rows[i]+' ['+snpApp.dt.y[i]+']')}
			// identify partition table
			snpApp.disp('computing partition table ...')
			snpApp.dt = snpCalc.partitions(snpApp.dt);
			snpApp.disp('ranked partitions:');
			I = jmat.array2mat(snpApp.dt.partition); // seprating binary indexes from ranksum values
			J=jmat.sort(I[1].map(function(ii){return -ii})); // sorting ranksums, descending
			snpApp.dt.Ind=J[1];
			snpApp.dt.partitionCodeList=I[0];
			snpApp.dt.rankSumList=I[1];
			var m = snpApp.dt.Ind.length;
			for(j=0;j<m;j++){ // look carefully here, it explain how to index ranked results
				var ii = snpApp.dt.Ind[j];
				snpApp.disp(j+'/'+m+' : '+snpApp.dt.partitionCodeList[ii]+' --> '+snpApp.dt.rankSumList[ii])
			}

		}
	},

	parseSnp:function(evt){
		//var linhas = evt.target.value.split('\n');
		var linhas = jmat.gId('SNPs').value.split('\n');
		snpApp.disp('SNP parsing ...');
		snpApp.dt={snps:[],rows:[],cols:[],y:[]};
		j=0;
		for (var i=0;i<linhas.length;i++){
			var linha = linhas[i].split('\t');
			snpApp.dt.snps[j]=linha[1];
			snpApp.dt.rows[j]=linha[0];
			//snpApp.dt.y[j]=linha[2].split(',').map(function(x){return JSON.parse(x)});
			j=j+1;	
		}
		var linhas = jmat.gId('phenotypes').value.split('\n');
		snpApp.disp('Phenotype parsing ...');
		j=0;
		for (var i=0;i<linhas.length;i++){
			snpApp.dt.y[j]=linhas[j].split(',').map(function(x){return JSON.parse(x)});
			j=j+1;	
		}
		
	}
};

snpCalc={

	partitions:function(dt){ // identifies p value for all possible partition profiles
		var n=dt.rows.length, pt=[];
		for(var i=1;i<=Math.pow(2,n)-2;i++){
			pt[i-1]=jmat.toBin(i,n);
		}
		dt.partition=[];
		for(var i=0;i<pt.length;i++){
			var x = [0]; // just to seed the concatenation
			var y = [0];
			for(var j=0;j<n;j++){
				if(pt[i][j]=='1'){x=jmat.cat(x,dt.y[j])}
				else{y=jmat.cat(y,dt.y[j])}
			}
			dt.partition[pt[i]]=jmat.ranksum(x.slice(1),y.slice(1)); //note seed sliced off
		}
		return dt
	}
}





snpApp.start();

