function UrlExists(url) {
    var http = new XMLHttpRequest();
    http.open('HEAD', url, false);
    http.send();
    return http.status!=404;
}

if (UrlExists('graphs/')) {
    $().ready(function(){
	var trArray = []
	$('#example tr').each(function(){
            var tdArray = []
            $(this).find('td').each(function(){
                tdArray.push($(this))
            })
		trArray.push(tdArray)     
	});
	var key;
	var guess;
	for(row = 0; row < trArray.length; row++){
	    // This will break things if the percents change :(
	    var interval = 0;
	    var thisrow = trArray[row];
	    for(cell = 8; cell < trArray[row].length; cell++){
		// Cell 0 is the peptide key
		key = trArray[row][0].text();	
		guess = trArray[row][cell].text();
		trArray[row][cell].html('<a href="graphs/' + key + '/' + interval + '.png" target="_blank">' + guess + '</a>');
		interval += 10;
	    }
	}
    })
} else {
    console.warn("graphs/ seems inaccessible. Not linking to graphs.");
}