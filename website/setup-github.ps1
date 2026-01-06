# PowerShell script to help set up GitHub repository
# Run this after creating the repository on GitHub

param(
    [Parameter(Mandatory=$true)]
    [string]$GitHubUsername,
    
    [Parameter(Mandatory=$false)]
    [string]$RepoName = "ignition-protocol-website"
)

Write-Host "Setting up GitHub repository..." -ForegroundColor Green
Write-Host ""

# Check if we're in the right directory
if (-not (Test-Path "package.json")) {
    Write-Host "Error: Please run this script from the website directory" -ForegroundColor Red
    exit 1
}

# Add remote
$remoteUrl = "https://github.com/$GitHubUsername/$RepoName.git"
Write-Host "Adding remote: $remoteUrl" -ForegroundColor Yellow
git remote add origin $remoteUrl

# Rename branch to main (if needed)
$currentBranch = git branch --show-current
if ($currentBranch -ne "main") {
    Write-Host "Renaming branch to 'main'..." -ForegroundColor Yellow
    git branch -M main
}

# Push to GitHub
Write-Host ""
Write-Host "Pushing to GitHub..." -ForegroundColor Yellow
Write-Host "You may be prompted for your GitHub credentials." -ForegroundColor Cyan
git push -u origin main

Write-Host ""
Write-Host "✅ Successfully pushed to GitHub!" -ForegroundColor Green
Write-Host ""
Write-Host "Next steps:" -ForegroundColor Cyan
Write-Host "1. Visit https://github.com/$GitHubUsername/$RepoName" -ForegroundColor White
Write-Host "2. Go to Vercel and import your repository for automatic deployment" -ForegroundColor White
Write-Host "3. Your website will be live at: https://$RepoName.vercel.app" -ForegroundColor White

