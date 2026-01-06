# Ignition Protocol Website

Website for the Ignition Protocol project—An Evolutionary Engine for Molecular Combustion.

## Overview

This is a Next.js 14 website that documents the fuel blend optimisation engine and showcases key highlights from the Ignition Protocol project book. The website features:

- **Engine Documentation**: Complete API reference, quick start guide, and usage examples
- **Project Highlights**: Curated excerpts from the project book
- **Interactive Demo**: Playground for optimising fuel blends
- **Flame & Darkness Theme**: Custom design language matching the "Ignition Protocol" aesthetic

## Tech Stack

- **Framework**: Next.js 14 (App Router)
- **Language**: TypeScript
- **Styling**: Tailwind CSS with custom flame/dark theme
- **Components**: shadcn/ui
- **Charts**: Recharts
- **Deployment**: Vercel

## Getting Started

### Development

```bash
# Install dependencies
npm install

# Run development server
npm run dev
```

Open [http://localhost:3000](http://localhost:3000) in your browser.

### Build

```bash
# Build for production
npm run build

# Start production server
npm start
```

## Project Structure

```
website/
├── app/                    # Next.js App Router pages
│   ├── engine/           # Engine documentation
│   ├── project/           # Project highlights
│   ├── demo/             # Interactive demo
│   └── api/              # API routes
├── components/            # React components
│   ├── ui/               # shadcn/ui components
│   ├── charts/           # Chart components
│   └── layout/           # Layout components
├── lib/                   # Utilities
├── public/                # Static assets
└── styles/               # Global styles
```

## Design System

The website utilises a custom "Flame & Darkness" theme:

- **Dark Foundation**: Deep black backgrounds (#0a0a0a, #1a1a1a)
- **Flame Accents**: Orange (#ff6b35), red-orange (#ff4500), yellow (#ffd700), amber (#ffb347)
- **Ember Glows**: Deep reds (#8b0000, #4a0000)
- **Animations**: Flame flicker, glow pulse, ember glow effects

## GitHub Repository Setup

The website is ready to be pushed to GitHub. Follow these steps:

### Quick Setup

1. **Create a new repository on GitHub:**
   - Navigate to https://github.com/new
   - Repository name: `ignition-protocol-website` (or your preferred name)
   - Description: "Website for Ignition Protocol—An Evolutionary Engine for Molecular Combustion"
   - Choose **Public** or **Private**
   - **DO NOT** initialise with README, .gitignore, or licence (these already exist)
   - Click "Create repository"

2. **Push your code:**
   ```powershell
   # Option 1: Utilise the setup script
   .\setup-github.ps1 -GitHubUsername YOUR_USERNAME
   
   # Option 2: Manual commands
   git remote add origin https://github.com/YOUR_USERNAME/ignition-protocol-website.git
   git branch -M main
   git push -u origin main
   ```

See `setup-github.md` for detailed instructions.

## Deployment

The website is configured for Vercel deployment. Simply connect your GitHub repository to Vercel and deploy.

### Environment Variables

No environment variables are required for basic functionality. The API routes may require Python engine access in production.

## Running on Linux (Old MacBook)

### Linux Distribution Recommendations

For running this website on an old MacBook with Linux:

**Pop!_OS** (Current choice—Recommended if working well):
- Excellent performance on older hardware
- Excellent developer tools support
- Straightforward to set up and maintain
- Good driver support for Mac hardware
- Based on Ubuntu, so excellent package availability

**Alternative Options:**

1. **Xubuntu** (If better performance is required):
   - XFCE desktop is more lightweight than PopOS
   - Better for very old hardware (2012–2015 MacBooks)
   - Still user-friendly and well-supported
   - Lower memory footprint

2. **Lubuntu** (Lightest option):
   - LXQt desktop—very lightweight
   - Best performance on old hardware
   - Minimal resource usage
   - Good for 2010–2013 MacBooks

3. **Fedora Workstation** (If cutting-edge features are desired):
   - Latest packages and features
   - Excellent for development
   - More resource-intensive than PopOS
   - Good for 2015+ MacBooks

**Recommendation**: 
- **If Pop!_OS is working well, maintain it**—It is an excellent choice for development and has good hardware support
- If performance issues are experienced, try **Xubuntu** (lighter desktop environment)
- For very old MacBooks (2010–2012), consider **Lubuntu** (lightest option)
- For newer MacBooks (2015+) requiring cutting-edge packages, **Fedora Workstation** is a solid alternative

**Why Pop!_OS is a good choice:**
- Based on Ubuntu LTS (long-term support, stable)
- Excellent NVIDIA GPU support (if your MacBook has one)
- Excellent for development work (Node.js, Python, etc.)
- Active community and good documentation
- System76's optimisations for performance

### Performance Tips for Old Hardware

1. **Utilise a lightweight desktop environment** (XFCE or LXQt)
2. **Disable unnecessary services** to free up RAM
3. **Utilise SSD** if possible (dramatically improves performance)
4. **Increase swap space** if limited RAM is available
5. **Close unnecessary applications** when running the dev server

## API Routes

- `/api/optimize` - Fuel blend optimisation endpoint
  - **Note**: Currently in demo mode—returns simulated results
  - For production, implement full Python engine integration
  - See `app/api/optimize/route.ts` for implementation details
- `/api/health` - Health check endpoint

## License

Developed for educational and research purposes.
