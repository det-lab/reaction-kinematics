document$.subscribe(() => {
    const table = $('#mass-table');

    if (table.length) {
        table.DataTable({
            pageLength: 25
        });
    }
});