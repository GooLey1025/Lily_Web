let blastnTable = null;
let blastpTable = null;

async function postForm(url, form) {
  const fd = new FormData(form);
  const resp = await fetch(url, { method: 'POST', body: fd });
  const data = await resp.json();
  if (!resp.ok) {
    throw new Error(data.error || 'Request failed');
  }
  return data;
}

function renderDataTable(tableId, rows) {
  const header = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'];
  
  if ($.fn.DataTable.isDataTable('#' + tableId)) {
    $('#' + tableId).DataTable().destroy();
  }
  
  $('#' + tableId).empty();
  
  const data = rows.map(r => header.map(h => r[h] ?? ''));
  
  const table = $('#' + tableId).DataTable({
    data: data,
    columns: header.map(h => ({ title: h })),
    pageLength: 25,
    dom: 'Bfrtip',
    buttons: [
      {
        extend: 'csv',
        text: 'Download filtered as CSV',
        filename: tableId + '_filtered',
        exportOptions: {
          modifier: {
            search: 'applied',
            order: 'applied'
          }
        }
      }
    ],
    order: [[10, 'asc']],
    language: {
      search: 'Filter:',
      lengthMenu: 'Show _MENU_ entries',
      info: 'Showing _START_ to _END_ of _TOTAL_ entries',
      infoFiltered: '(filtered from _MAX_ total entries)'
    }
  });
  
  return table;
}

function showError(errorId, message) {
  const errorEl = document.getElementById(errorId);
  errorEl.textContent = message;
  errorEl.style.display = 'block';
}

function hideError(errorId) {
  const errorEl = document.getElementById(errorId);
  errorEl.style.display = 'none';
}

function showLoading(loadingId) {
  const loadingEl = document.getElementById(loadingId);
  loadingEl.style.display = 'block';
}

function hideLoading(loadingId) {
  const loadingEl = document.getElementById(loadingId);
  loadingEl.style.display = 'none';
}

document.getElementById('form-blastn').addEventListener('submit', async (e) => {
  e.preventDefault();
  hideError('error-blastn');
  showLoading('loading-blastn');
  try {
    const data = await postForm('/blastn', e.target);
    blastnTable = renderDataTable('table-blastn', data.rows);
    const dl = document.getElementById('dl-blastn');
    dl.href = data.download.results;
    dl.style.display = 'inline-block';
  } catch (err) {
    showError('error-blastn', 'Error: ' + err.message);
  } finally {
    hideLoading('loading-blastn');
  }
});

document.getElementById('form-blastp').addEventListener('submit', async (e) => {
  e.preventDefault();
  hideError('error-blastp');
  showLoading('loading-blastp');
  try {
    const data = await postForm('/blastp', e.target);
    blastpTable = renderDataTable('table-blastp', data.rows);
    const dl = document.getElementById('dl-blastp');
    dl.href = data.download.results;
    dl.style.display = 'inline-block';
  } catch (err) {
    showError('error-blastp', 'Error: ' + err.message);
  } finally {
    hideLoading('loading-blastp');
  }
});

document.getElementById('form-subseq').addEventListener('submit', async (e) => {
  e.preventDefault();
  hideError('error-subseq');
  showLoading('loading-subseq');
  try {
    const data = await postForm('/subseq', e.target);
    const dl = document.getElementById('dl-subseq');
    dl.href = data.download.fasta;
    dl.style.display = 'inline-block';
  } catch (err) {
    showError('error-subseq', 'Error: ' + err.message);
  } finally {
    hideLoading('loading-subseq');
  }
});


