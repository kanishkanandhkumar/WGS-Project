/**
 * OmniVariant Analytics Suite - Application Engine
 */

let vcfFiles = [];
let masterVariants = [];
let activeCharts = {};
let detectedPipeline = "None Detected";

// DOM Elements
const vcfInput = document.getElementById('vcf-files');
const dropZone = document.getElementById('drop-zone');
const loadBtn = document.getElementById('load-btn');
const filterGene = document.getElementById('filter-gene');
const filterConsequence = document.getElementById('filter-consequence');
const filterImpact = document.getElementById('filter-impact');
const filterChrom = document.getElementById('filter-chrom');
const exportCsvBtn = document.getElementById('export-csv-btn');
const resetBtn = document.getElementById('reset-btn');
const modal = document.getElementById('variant-modal');
const closeModalBtn = document.getElementById('close-modal');

// --- FILE DRAG & DROP HANDLERS ---
dropZone.addEventListener('dragover', (e) => {
    e.preventDefault();
    dropZone.classList.add('dropzone-active');
});
dropZone.addEventListener('dragleave', () => dropZone.classList.remove('dropzone-active'));
dropZone.addEventListener('drop', (e) => {
    e.preventDefault();
    dropZone.classList.remove('dropzone-active');
    if (e.dataTransfer.files.length) {
        vcfInput.files = e.dataTransfer.files;
        handleFileSelection(Array.from(e.dataTransfer.files));
    }
});

vcfInput.addEventListener('change', (e) => handleFileSelection(Array.from(e.target.files)));

function handleFileSelection(files) {
    vcfFiles = files;
    const statusEl = document.getElementById('vcf-status');
    statusEl.classList.remove('hidden');
    statusEl.innerText = `📦 ${vcfFiles.length} file(s) attached.`;
}

resetBtn.addEventListener('click', () => location.reload());

// --- REAL-TIME FILTER ENGINE ---
const applyFilters = () => {
    const geneQuery = filterGene.value.toLowerCase().trim();
    const consequenceQuery = filterConsequence.value.toLowerCase().trim();
    const impactQuery = filterImpact.value;
    const chromQuery = filterChrom.value;

    const filteredData = masterVariants.filter(v => {
        const matchesGene = v.gene.toLowerCase().includes(geneQuery);
        const matchesConsequence = v.consequence.toLowerCase().includes(consequenceQuery);
        const matchesImpact = impactQuery === 'ALL' || v.impactLevel === impactQuery;
        const matchesChrom = chromQuery === 'ALL' || v.chrom === chromQuery;
        return matchesGene && matchesConsequence && matchesImpact && matchesChrom;
    });

    renderTable(filteredData.slice(0, 500), filteredData.length);
};

[filterGene, filterConsequence, filterImpact, filterChrom].forEach(el => el.addEventListener('input', applyFilters));

// --- PARSING CONTROLLER ---
loadBtn.addEventListener('click', async () => {
    if (vcfFiles.length === 0) { 
        alert("Please select or drop a VCF file first."); 
        return; 
    }
    
    document.getElementById('welcome-state').classList.add('hidden');
    document.getElementById('loading-state').classList.remove('hidden');
    await new Promise(resolve => setTimeout(resolve, 50));

    try {
        const fastMode = document.getElementById('fast-mode').checked;
        masterVariants = []; 
        let pipelineScores = { CSQ: 0, ANN: 0 };
        
        let stats = { 
            total: 0, snp: 0, mapped: 0, high_impact: 0, transitions: 0, transversions: 0,
            chr: {}, 
            types: { SNV: 0, INS: 0, DEL: 0, MNP: 0, COMPLEX: 0 },
            impact: { HIGH: 0, MODERATE: 0, LOW: 0, MODIFIER: 0 }
        };

        for (let file of vcfFiles) {
            const fileContent = await readFileContent(file);
            parseVCFStream(fileContent, stats, fastMode, pipelineScores);
        }

        detectedPipeline = pipelineScores.CSQ >= pipelineScores.ANN && pipelineScores.CSQ > 0 
            ? "Ensembl VEP (CSQ)" 
            : pipelineScores.ANN > 0 ? "SnpEff (ANN)" : "Raw Native VCF";

        populateChromosomeDropdown(Object.keys(stats.chr));
        renderDashboard(stats);
    } catch (err) {
        console.error(err);
        alert("An error occurred while parsing the VCF file.");
        document.getElementById('loading-state').classList.add('hidden');
        document.getElementById('welcome-state').classList.remove('hidden');
    }
});

function readFileContent(file) {
    return new Promise((resolve, reject) => { 
        const reader = new FileReader(); 
        reader.onload = e => resolve(e.target.result); 
        reader.onerror = err => reject(err); 
        reader.readAsText(file); 
    });
}

// --- OPTIMIZED VCF STREAM PARSER ---
function parseVCFStream(vcfText, stats, fastMode, pipelineScores) {
    const loopLimit = fastMode ? 25000 : 1000000; 
    let parsedLines = 0;
    let currentPos = 0;
    let nextNL = 0;

    while ((nextNL = vcfText.indexOf('\n', currentPos)) !== -1 && parsedLines < loopLimit) {
        const rawLine = vcfText.substring(currentPos, nextNL).trim();
        currentPos = nextNL + 1;

        if (!rawLine || rawLine.startsWith('#')) continue;

        stats.total++;
        parsedLines++;
        
        const cols = rawLine.split('\t');
        if (cols.length < 8) continue;

        const chrom = cols[0];
        const pos = parseInt(cols[1]);
        const ref = cols[3].toUpperCase();
        const alt = cols[4].toUpperCase();
        const info = cols[7];

        stats.chr[chrom] = (stats.chr[chrom] || 0) + 1;

        // Classify Variant Type & Ti/Tv
        let varType = "COMPLEX";
        if (ref.length === 1 && alt.length === 1) {
            varType = "SNV";
            stats.snp++;
            if (isTransition(ref, alt)) stats.transitions++;
            else stats.transversions++;
        } else if (ref.length < alt.length) {
            varType = "INS";
        } else if (ref.length > alt.length) {
            varType = "DEL";
        } else if (ref.length === alt.length) {
            varType = "MNP";
        }
        stats.types[varType] = (stats.types[varType] || 0) + 1;

        let gene = "-";
        let consequence = "-";
        let impactLevel = "MODIFIER";

        if (info.includes('CSQ=')) {
            pipelineScores.CSQ++;
            const top = info.split('CSQ=')[1].split(';')[0].split(',')[0];
            const fields = top.split('|');
            if (fields.length > 3) {
                consequence = fields[1] ? fields[1].replace(/_/g, ' ') : "unknown";
                impactLevel = fields[2] || "MODIFIER";
                gene = fields[3] || "-";
            }
        } else if (info.includes('ANN=')) {
            pipelineScores.ANN++;
            const top = info.split('ANN=')[1].split(';')[0].split(',')[0];
            const fields = top.split('|');
            if (fields.length > 3) {
                consequence = fields[1] ? fields[1].replace(/_/g, ' ') : "unknown";
                impactLevel = fields[2] || "MODIFIER";
                gene = fields[3] || "-";
            }
        }

        if (gene !== "-") stats.mapped++;
        
        const normImpact = impactLevel.toUpperCase().trim();
        if (stats.impact[normImpact] !== undefined) {
            stats.impact[normImpact]++;
            if (normImpact === 'HIGH') stats.high_impact++;
        } else {
            stats.impact['MODIFIER']++;
        }

        masterVariants.push({ chrom, pos, ref, alt, varType, gene, consequence, impactLevel: normImpact, rawInfo: info });
    }
}

function isTransition(ref, alt) {
    const pair = ref + alt;
    return pair === "AG" || pair === "GA" || pair === "CT" || pair === "TC";
}

function populateChromosomeDropdown(contigs) {
    filterChrom.innerHTML = '<option value="ALL">ALL CONTIGS</option>';
    contigs.sort().forEach(c => {
        const opt = document.createElement('option');
        opt.value = c;
        opt.innerText = c;
        filterChrom.appendChild(opt);
    });
}

// --- DASHBOARD RENDERER ---
function renderDashboard(stats) {
    document.getElementById('metric-total').innerText = stats.total.toLocaleString();
    document.getElementById('metric-snps').innerText = stats.snp.toLocaleString();
    document.getElementById('metric-high').innerText = stats.high_impact.toLocaleString();
    document.getElementById('metric-mapped').innerText = stats.mapped.toLocaleString();

    const tiTv = stats.transversions > 0 ? (stats.transitions / stats.transversions).toFixed(2) : "N/A";
    document.getElementById('metric-titv').innerText = tiTv;

    const badge = document.getElementById('pipeline-badge');
    badge.innerText = detectedPipeline;
    badge.className = detectedPipeline.includes('VEP') 
        ? "text-[9px] font-mono font-bold tracking-widest uppercase px-2 py-0.5 rounded border border-indigo-500/30 text-indigo-400 bg-indigo-500/5"
        : detectedPipeline.includes('SnpEff')
        ? "text-[9px] font-mono font-bold tracking-widest uppercase px-2 py-0.5 rounded border border-teal-500/30 text-teal-400 bg-teal-500/5"
        : "text-[9px] font-mono font-bold tracking-widest uppercase px-2 py-0.5 rounded border border-slate-700 text-slate-400 bg-slate-800";

    [filterGene, filterConsequence, filterImpact, filterChrom].forEach(el => el.disabled = false);

    renderCharts(stats);
    renderTable(masterVariants.slice(0, 500), masterVariants.length);

    document.getElementById('loading-state').classList.add('hidden');
    document.getElementById('dashboard-view').classList.remove('hidden');
}

// --- CHARTS GENERATOR ---
function renderCharts(stats) {
    // 1. Contig Bar Chart
    if (activeCharts.chr) activeCharts.chr.destroy();
    const sortedContigs = Object.entries(stats.chr).sort((a,b) => b[1] - a[1]).slice(0, 15);
    activeCharts.chr = new Chart(document.getElementById('chrChart'), {
        type: 'bar',
        data: {
            labels: sortedContigs.map(i => i[0]),
            datasets: [{ data: sortedContigs.map(i => i[1]), backgroundColor: '#14b8a6', borderRadius: 4 }]
        },
        options: {
            responsive: true, maintainAspectRatio: false,
            plugins: { legend: { display: false } },
            scales: {
                x: { grid: { display: false }, ticks: { color: '#64748b', font: { family: 'monospace', size: 9 } } },
                y: { grid: { color: '#1e293b' }, ticks: { color: '#64748b', font: { family: 'monospace', size: 9 } } }
            }
        }
    });

    // 2. Impact Doughnut Chart
    if (activeCharts.impact) activeCharts.impact.destroy();
    activeCharts.impact = new Chart(document.getElementById('impactChart'), {
        type: 'doughnut',
        data: {
            labels: Object.keys(stats.impact),
            datasets: [{ data: Object.values(stats.impact), backgroundColor: ['#f43f5e', '#f59e0b', '#3b82f6', '#475569'], borderWidth: 0 }]
        },
        options: {
            responsive: true, maintainAspectRatio: false, cutout: '75%',
            plugins: { legend: { position: 'right', labels: { color: '#94a3b8', font: { family: 'monospace', size: 9 }, boxWidth: 10 } } }
        }
    });

    // 3. Variant Type Pie Chart
    if (activeCharts.type) activeCharts.type.destroy();
    activeCharts.type = new Chart(document.getElementById('typeChart'), {
        type: 'pie',
        data: {
            labels: Object.keys(stats.types),
            datasets: [{ data: Object.values(stats.types), backgroundColor: ['#10b981', '#6366f1', '#ec4899', '#8b5cf6', '#64748b'], borderWidth: 0 }]
        },
        options: {
            responsive: true, maintainAspectRatio: false,
            plugins: { legend: { position: 'right', labels: { color: '#94a3b8', font: { family: 'monospace', size: 9 }, boxWidth: 10 } } }
        }
    });
}

// --- MATRIX TABLE & MODAL INSPECTOR ---
function renderTable(rows, totalCount) {
    const tbody = document.getElementById('variant-rows');
    document.getElementById('table-counter').innerText = `Showing up to 500 of ${totalCount.toLocaleString()} records`;
    
    if (rows.length === 0) { 
        tbody.innerHTML = `<tr><td colspan="8" class="p-6 text-center text-slate-600 italic">No matching variant annotations found.</td></tr>`; 
        return; 
    }
    
    tbody.innerHTML = rows.map((v, idx) => {
        let badgeClass = 'badge-modifier';
        if (v.impactLevel === 'HIGH') badgeClass = 'badge-high';
        else if (v.impactLevel === 'MODERATE') badgeClass = 'badge-moderate';
        else if (v.impactLevel === 'LOW') badgeClass = 'badge-low';

        return `
        <tr class="hover:bg-slate-900/60 cursor-pointer transition-colors" onclick="openVariantModal(${idx})">
            <td class="p-3 text-slate-400 font-semibold">${v.chrom}</td>
            <td class="p-3 text-slate-300 font-medium">${v.pos.toLocaleString()}</td>
            <td class="p-3 text-teal-400 font-bold">${v.ref}</td>
            <td class="p-3 text-indigo-400 font-bold">${v.alt}</td>
            <td class="p-3 text-slate-400 font-mono text-[10px]">${v.varType}</td>
            <td class="p-3 text-slate-200 font-bold tracking-wide">${v.gene}</td>
            <td class="p-3 text-slate-400 lowercase truncate max-w-xs">${v.consequence}</td>
            <td class="p-3 text-center"><span class="inline-block px-2 py-0.5 rounded text-[9px] font-bold ${badgeClass}">${v.impactLevel}</span></td>
        </tr>
    `}).join('');
}

window.openVariantModal = (index) => {
    const v = masterVariants[index];
    if (!v) return;

    document.getElementById('modal-title').innerText = `${v.chrom}:${v.pos} (${v.ref} → ${v.alt})`;
    document.getElementById('modal-body').innerHTML = `
        <p><span class="text-slate-500">Gene Symbol:</span> <strong class="text-emerald-400">${v.gene}</strong></p>
        <p><span class="text-slate-500">Consequence:</span> ${v.consequence}</p>
        <p><span class="text-slate-500">Impact Level:</span> ${v.impactLevel}</p>
        <p><span class="text-slate-500">Variant Class:</span> ${v.varType}</p>
        <div class="mt-3 pt-3 border-t border-slate-800">
            <span class="text-slate-500 block mb-1">INFO Field Snippet:</span>
            <div class="bg-slate-950 p-2 rounded text-[10px] break-all max-h-32 overflow-y-auto text-slate-400 font-mono">${v.rawInfo}</div>
        </div>
    `;
    modal.classList.remove('hidden');
};

closeModalBtn.addEventListener('click', () => modal.classList.add('hidden'));

// --- CSV EXPORT FUNCTION ---
exportCsvBtn.addEventListener('click', () => {
    if (!masterVariants.length) return;
    
    const headers = ["Contig", "Position", "Ref", "Alt", "Type", "Gene", "Consequence", "Impact"];
    const rows = masterVariants.map(v => [v.chrom, v.pos, v.ref, v.alt, v.varType, v.gene, v.consequence, v.impactLevel]);
    
    const csvContent = "data:text/csv;charset=utf-8," + [headers.join(","), ...rows.map(e => e.join(","))].join("\n");
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", "omnivariant_export.csv");
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
});
