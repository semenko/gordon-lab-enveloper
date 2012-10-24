/* Table initialisation */
$(document).ready(function() {
    oTable = $('#example').dataTable({
	/* Edge to edge without:
	"sDom": "<'row'<'span6'l><'span6'f>r>t<'row'<'span6'i><'span6'p>>",
	*/
	"sPaginationType": "bootstrap",
	"oLanguage": {
	    "sLengthMenu": "_MENU_ records per page"
	},
        "fnDrawCallback": function ( oSettings ) {
            if ( oSettings.aiDisplay.length == 0 )
            {
                return;
            }
             
            var nTrs = $('#example tbody tr');
            var iColspan = nTrs[0].getElementsByTagName('td').length;
            var sLastGroup = "";
            for ( var i=0 ; i<nTrs.length ; i++ )
            {
                var iDisplayIndex = oSettings._iDisplayStart + i;
                var sGroup = oSettings.aoData[ oSettings.aiDisplay[iDisplayIndex] ]._aData[0];
                if ( sGroup != sLastGroup )
                {
                    var nGroup = document.createElement( 'tr' );
                    var nCell = document.createElement( 'td' );
                    nCell.colSpan = iColspan;
                    nCell.className = "group";
                    nCell.innerHTML = sGroup;
                    nGroup.appendChild( nCell );
                    nTrs[i].parentNode.insertBefore( nGroup, nTrs[i] );
                    sLastGroup = sGroup;
                }
            }
        },
        "aoColumnDefs": [
            { "bVisible": false, "aTargets": [ 0 ] }
        ],
        "aaSortingFixed": [[ 0, 'asc' ]],
        "aaSorting": [[ 1, 'asc' ]],
        "sDom": 'lfr<"giveHeight"t>ip',
	"bStateSave": true,
        "fnStateSave": function (oSettings, oData) {
            localStorage.setItem('dt-' + window.location.pathname, JSON.stringify(oData));
        },
        "fnStateLoad": function (oSettings) {
            return JSON.parse(localStorage.getItem('dt-' + window.location.pathname));
        }
    });
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