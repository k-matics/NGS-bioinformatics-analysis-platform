/* Front-end helpers. No framework: plain JS + Bootstrap. All server text is
   HTML-escaped before insertion. */
const NGS = (() => {
  const BADGE = { pending: 'secondary', queued: 'info', running: 'primary', paused: 'warning',
                  completed: 'success', failed: 'danger', cancelled: 'dark', skipped: 'light' };

  const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) =>
    ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

  const badge = (status) =>
    `<span class="badge text-bg-${BADGE[status] || 'secondary'} status-badge">${esc(status).toUpperCase()}</span>`;

  const progressBar = (pct, status) =>
    `<div class="progress" style="height:1.1rem"><div class="progress-bar ${status === 'running' ? 'progress-bar-striped progress-bar-animated' : ''} bg-${BADGE[status] || 'primary'}" style="width:${pct}%">${pct}%</div></div>`;

  const duration = (sec) => {
    if (sec === null || sec === undefined || sec === '') return '—';
    sec = Math.floor(Number(sec));
    const h = Math.floor(sec / 3600), m = Math.floor((sec % 3600) / 60), s = sec % 60;
    const mm = String(m).padStart(2, '0'), ss = String(s).padStart(2, '0');
    return h ? `${h}h ${mm}m ${ss}s` : `${m}m ${ss}s`;
  };

  const meter = (label, pct, detail) => {
    const color = pct > 90 ? 'danger' : pct > 75 ? 'warning' : 'success';
    return `<div class="resource-meter mb-3"><div class="d-flex justify-content-between small"><span>${esc(label)}</span><span>${esc(detail)}</span></div>
      <div class="progress"><div class="progress-bar bg-${color}" style="width:${pct}%"></div></div></div>`;
  };

  async function pollSystem(el) {
    const load = async () => {
      try {
        const d = await (await fetch('/api/system')).json();
        let html = meter('CPU', d.cpu_percent, `${d.cpu_percent}% of ${d.cpu_count} threads · load ${d.load_avg.join(' / ')}`);
        html += meter('Memory', d.memory.used_percent, `${d.memory.available_gb} GB free of ${d.memory.total_gb} GB`);
        for (const disk of d.disks) {
          if (disk.same_device_as_previous) { html += `<div class="small text-secondary mb-2">${esc(disk.label)}: same disk as above</div>`; continue; }
          html += meter(disk.label, disk.used_percent, `${disk.free_gb} GB free of ${disk.total_gb} GB`);
        }
        el.innerHTML = html;
      } catch (e) { el.innerHTML = '<span class="text-danger small">Could not load server stats.</span>'; }
    };
    load(); setInterval(load, 5000);
  }

  async function pollDashboard(tbody) {
    const load = async () => {
      try {
        const d = await (await fetch('/api/dashboard')).json();
        const c = d.counts;
        const set = (k, v) => { const n = document.querySelector(`[data-count="${k}"]`); if (n) n.textContent = v; };
        set('queued', c.pending + c.queued); set('running', c.running + c.paused);
        set('completed', c.completed); set('failed', c.failed); set('cancelled', c.cancelled);
        tbody.innerHTML = d.active.length ? d.active.map((r) => `<tr>
            <td><a href="${esc(r.url)}">${esc(r.run_id)}</a></td><td>${esc(r.sample)}</td><td>${esc(r.workflow)}</td>
            <td>${badge(r.status)}</td><td>${esc(r.current_step || '—')}</td><td>${progressBar(r.progress, r.status)}</td></tr>`).join('')
          : '<tr><td colspan="6" class="text-secondary text-center py-3">No active jobs.</td></tr>';
      } catch (e) { /* keep last view */ }
    };
    setInterval(load, 5000);
  }

  function fileBrowser(modalId, inputId, startPath) {
    const modal = document.getElementById(modalId);
    const list = document.getElementById('browser-list');
    const pathEl = document.getElementById('browser-path');
    const input = document.getElementById(inputId);
    const open = async (path) => {
      const r = await fetch('/api/browse?path=' + encodeURIComponent(path));
      const d = await r.json();
      if (d.error) { list.innerHTML = `<li class="list-group-item text-danger">${esc(d.error)}</li>`; return; }
      pathEl.textContent = d.path;
      let html = d.parent ? `<li class="list-group-item list-group-item-action" data-dir="${esc(d.parent)}"><i class="bi bi-arrow-up"></i> ..</li>` : '';
      html += d.entries.map((e) => e.is_dir
        ? `<li class="list-group-item list-group-item-action" data-dir="${esc(e.path)}"><i class="bi bi-folder"></i> ${esc(e.name)}</li>`
        : `<li class="list-group-item list-group-item-action d-flex justify-content-between" data-file="${esc(e.path)}"><span><i class="bi bi-file-earmark"></i> ${esc(e.name)}</span><span class="text-secondary small">${(e.size / 1024 ** 3).toFixed(2)} GB</span></li>`).join('');
      list.innerHTML = html || '<li class="list-group-item text-secondary">No NGS files here.</li>';
    };
    list.addEventListener('click', (ev) => {
      const li = ev.target.closest('li'); if (!li) return;
      if (li.dataset.dir) open(li.dataset.dir);
      if (li.dataset.file) { input.value = li.dataset.file; bootstrap.Modal.getInstance(modal).hide(); }
    });
    modal.addEventListener('show.bs.modal', () => open(input.value || startPath));
  }

  function monitorRun(opts) {
    const logs = {};           // step key -> text
    let selected = null;       // null = all steps
    let elapsed = null, status = opts.initialStatus;
    const panel = document.getElementById('log-panel');
    const follow = document.getElementById('log-follow');

    const renderLog = () => {
      panel.textContent = selected ? (logs[selected] || '(no log yet)') : Object.values(logs).join('\n');
      if (follow.checked) panel.scrollTop = panel.scrollHeight;
    };
    document.querySelectorAll('.step-row').forEach((row) => row.addEventListener('click', () => {
      const key = row.dataset.step;
      selected = selected === key ? null : key;
      document.querySelectorAll('.step-row').forEach((r) => r.classList.toggle('selected', r.dataset.step === selected));
      document.getElementById('log-title').textContent = selected ? row.cells[1].firstChild.textContent.trim() : 'all steps';
      renderLog();
    }));

    setInterval(() => {
      if (elapsed !== null && (status === 'running' || status === 'paused')) {
        elapsed += 1; document.getElementById('run-elapsed').textContent = duration(elapsed);
      }
    }, 1000);

    const es = new EventSource(opts.eventsUrl);
    es.addEventListener('state', (ev) => {
      const s = JSON.parse(ev.data);
      status = s.status; elapsed = s.elapsed;
      document.getElementById('run-status').innerHTML = badge(s.status);
      document.getElementById('run-progress').innerHTML = progressBar(s.progress, s.status);
      document.getElementById('run-progress-label').textContent = s.progress + '%';
      document.getElementById('run-elapsed').textContent = duration(s.elapsed);
      document.getElementById('current-step').textContent = s.current_step || '—';
      document.getElementById('current-command').textContent = s.current_command || '—';
      const err = document.getElementById('run-error');
      err.textContent = s.error || ''; err.classList.toggle('d-none', !s.error);
      const warn = document.getElementById('run-warnings');
      if (s.warnings && s.warnings.length) {
        warn.innerHTML = '<strong>Warnings (need review, not necessarily errors):</strong><ul class="mb-0">' +
          s.warnings.map((w) => `<li>${esc(w)}</li>`).join('') + '</ul>';
        warn.classList.remove('d-none');
      }
      for (const st of s.steps) {
        const row = document.querySelector(`.step-row[data-step="${st.key}"]`);
        if (!row) continue;
        row.querySelector('.step-status').innerHTML = badge(st.status) + (st.warnings ? ` <span class="badge text-bg-warning">${st.warnings} warn</span>` : '');
        row.querySelector('.step-duration').textContent = duration(st.duration);
        row.querySelector('.step-msg').textContent = st.message || '';
      }
      document.querySelectorAll('[data-show-when]').forEach((f) => {
        f.hidden = !f.dataset.showWhen.split(',').includes(s.status);
      });
    });
    es.addEventListener('log', (ev) => {
      const d = JSON.parse(ev.data);
      logs[d.step] = (logs[d.step] || '') + d.text;
      if (!selected || selected === d.step) renderLog();
    });
    es.addEventListener('end', () => {
      es.close();
      // Re-render once so outputs, commands and audit trail are shown in full.
      if (!opts.terminal.includes(opts.initialStatus)) setTimeout(() => location.reload(), 800);
    });
  }

  // Start-workflow form: show the selected workflow's parameters, pre-select
  // matching input files, show per-caller parameters only for ticked callers.
  function workflowForm(selectId) {
    const sel = document.getElementById(selectId);
    if (!sel) return;
    const form = sel.form;
    function syncCallers() {
      form.querySelectorAll('input[data-toggle-target]').forEach(cb => {
        const box = document.getElementById(cb.dataset.toggleTarget);
        if (box) box.style.display = cb.checked ? '' : 'none';
      });
    }
    function update() {
      const wf = sel.value;
      form.querySelectorAll('.wf-panel').forEach(p => {
        const on = p.dataset.wf === wf;
        p.style.display = on ? '' : 'none';
        p.querySelectorAll('input,select').forEach(i => { i.disabled = !on || i.dataset.unavailable === '1'; });
      });
      form.querySelectorAll('.wf-desc').forEach(d => { d.style.display = d.dataset.wf === wf ? '' : 'none'; });
      const roles = (sel.selectedOptions[0].dataset.roles || '').split(',').filter(Boolean);
      const checks = form.querySelectorAll('.file-check');
      let bams = 0;
      checks.forEach(c => {
        let on = roles.length === 0 || roles.includes(c.dataset.role);
        if (on && c.dataset.role === 'bam' && roles.length) { bams += 1; on = bams === 1; }  // one BAM per BAM->VCF run
        c.checked = on;
      });
      syncCallers();
    }
    form.querySelectorAll('input[data-toggle-target]').forEach(cb => {
      if (cb.disabled) cb.dataset.unavailable = '1';
      cb.addEventListener('change', syncCallers);
    });
    sel.addEventListener('change', update);
    update();
  }


  // Parameter presets on the Start form: fill the fields of the selected workflow from a preset.
  // Only the fields a preset contains are changed; the user can still edit everything afterwards.
  function presetPicker(selectId, wfSelectId, presets, infoId, hiddenId) {
    const ps = document.getElementById(selectId), wfSel = document.getElementById(wfSelectId);
    if (!ps || !wfSel) return;
    const form = wfSel.form, info = document.getElementById(infoId), hidden = document.getElementById(hiddenId);
    const byId = {}; presets.forEach(p => { byId[p.id] = p; });
    function setField(name, value) {
      const els = form.querySelectorAll(`[name="${CSS.escape(name)}"]`);
      if (!els.length) return false;
      els.forEach(el => {
        if (el.type === 'checkbox' && Array.isArray(value)) {
          el.checked = el.dataset.unavailable !== '1' && value.map(String).includes(el.value);
        } else if (el.type === 'checkbox') { el.checked = !!value; }
        else { el.value = value; }
        el.dispatchEvent(new Event('change'));
      });
      return true;
    }
    function apply() {
      const p = byId[ps.value];
      hidden.value = p ? p.id : '';
      if (!p) { info.textContent = ''; return; }
      const wf = wfSel.value; let set = 0, skipped = 0;
      Object.entries(p.values).forEach(([k, v]) => {
        const name = k.startsWith('p:') ? `p.${wf}.${k.slice(2)}` : `c.${wf}.${k.slice(2)}`;
        if (setField(name, v)) set += 1; else skipped += 1;
      });
      info.textContent = `Preset "${p.name}" applied: ${set} field(s) set` +
        (skipped ? `, ${skipped} not used by this workflow` : '') + '. You can still change any value; the run records what you changed.';
    }
    ps.addEventListener('change', apply);
    wfSel.addEventListener('change', () => { if (ps.value) apply(); });
  }

  return { pollSystem, pollDashboard, fileBrowser, monitorRun, workflowForm, presetPicker, esc };
})();
