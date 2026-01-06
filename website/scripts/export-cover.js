const puppeteer = require('puppeteer');
const path = require('path');
const fs = require('fs');

async function exportBookCover() {
  console.log('Starting book cover export...');
  
  const browser = await puppeteer.launch({
    headless: true,
    args: ['--no-sandbox', '--disable-setuid-sandbox']
  });
  
  try {
    const page = await browser.newPage();
    
    // Set viewport to match book cover size (8.5" × 11" at 300 DPI)
    // 8.5 inches × 300 DPI = 2550 pixels width
    // 11 inches × 300 DPI = 3300 pixels height
    await page.setViewport({
      width: 2550,
      height: 3300,
      deviceScaleFactor: 1
    });
    
    // Get the absolute path to the HTML file
    const htmlPath = path.join(__dirname, '../public/book-cover.html');
    const fileUrl = `file://${htmlPath}`;
    
    console.log(`Loading: ${fileUrl}`);
    await page.goto(fileUrl, {
      waitUntil: 'networkidle0',
      timeout: 30000
    });
    
    // Wait a bit for any animations/rendering
    await page.waitForTimeout(1000);
    
    // Create output directory if it doesn't exist
    const outputDir = path.join(__dirname, '../public/exports');
    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
    }
    
    // Export as PNG
    const outputPath = path.join(outputDir, 'book-cover.png');
    await page.screenshot({
      path: outputPath,
      type: 'png',
      fullPage: true,
      printBackground: true
    });
    
    console.log(`✅ Book cover exported to: ${outputPath}`);
    console.log(`   Size: 2550 × 3300 pixels (8.5" × 11" at 300 DPI)`);
    
  } catch (error) {
    console.error('Error exporting book cover:', error);
    process.exit(1);
  } finally {
    await browser.close();
  }
}

// Also export poster
async function exportPoster() {
  console.log('\nStarting poster export...');
  
  const browser = await puppeteer.launch({
    headless: true,
    args: ['--no-sandbox', '--disable-setuid-sandbox']
  });
  
  try {
    const page = await browser.newPage();
    
    // Set viewport to match poster size (24" × 36" at 300 DPI)
    // 24 inches × 300 DPI = 7200 pixels width
    // 36 inches × 300 DPI = 10800 pixels height
    await page.setViewport({
      width: 7200,
      height: 10800,
      deviceScaleFactor: 1
    });
    
    const htmlPath = path.join(__dirname, '../public/poster.html');
    const fileUrl = `file://${htmlPath}`;
    
    console.log(`Loading: ${fileUrl}`);
    await page.goto(fileUrl, {
      waitUntil: 'networkidle0',
      timeout: 30000
    });
    
    await page.waitForTimeout(1000);
    
    const outputDir = path.join(__dirname, '../public/exports');
    if (!fs.existsSync(outputDir)) {
      fs.mkdirSync(outputDir, { recursive: true });
    }
    
    const outputPath = path.join(outputDir, 'poster.png');
    await page.screenshot({
      path: outputPath,
      type: 'png',
      fullPage: true,
      printBackground: true
    });
    
    console.log(`✅ Poster exported to: ${outputPath}`);
    console.log(`   Size: 7200 × 10800 pixels (24" × 36" at 300 DPI)`);
    
  } catch (error) {
    console.error('Error exporting poster:', error);
    process.exit(1);
  } finally {
    await browser.close();
  }
}

// Run exports
(async () => {
  await exportBookCover();
  await exportPoster();
  console.log('\n✨ All exports complete!');
})();














