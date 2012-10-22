/* Table initialisation */
$(document).ready(function() {
    $('#example').dataTable( {
	/* Edge to edge without:
	"sDom": "<'row'<'span6'l><'span6'f>r>t<'row'<'span6'i><'span6'p>>",
	*/
	"sPaginationType": "bootstrap",
	"oLanguage": {
	    "sLengthMenu": "_MENU_ records per page"
	}
    } );
    // Anything starting with a #
    var anchor = document.URL.split('#')[1];
    if (anchor) {
	$('#example').dataTable().fnFilter(anchor);
	// Set the search form, too. XSS here?
	$('#example_filter input').val(anchor);
    };
} );
