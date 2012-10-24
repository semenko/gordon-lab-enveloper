/* Table initialisation */
$(document).ready(function() {
    $('#example').dataTable( {
	/* Edge to edge without:
	"sDom": "<'row'<'span6'l><'span6'f>r>t<'row'<'span6'i><'span6'p>>",
	*/
	"sPaginationType": "bootstrap",
	"oLanguage": {
	    "sLengthMenu": "_MENU_ records per page"
	},
	"bStateSave": true,
        "fnStateSave": function (oSettings, oData) {
            localStorage.setItem('dt-' + window.location.pathname, JSON.stringify(oData));
        },
        "fnStateLoad": function (oSettings) {
            return JSON.parse(localStorage.getItem('dt-' + window.location.pathname));
        }
    } );
    // Anything starting with a #
    function updateHashFilter() {
	var anchor = document.URL.split('#')[1];
	if (anchor) {
	    $('#example').dataTable().fnFilter(anchor);
	    // Set the search form, too. XSS here?
	    $('#example_filter input').val(anchor);
	};
    };
    updateHashFilter();
    $(window).bind('hashchange', function() {
	updateHashFilter();
    });
} );
